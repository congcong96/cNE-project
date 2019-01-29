function nespktrain = ne_split_spon_evoked_spkstruct(spk, trigger, stimlen)

% Splits spon-evoked sorted trains. Assumes that trigger is only for evoked
% activity, and will consider all data up to that point to be spontaneous
% activity.

rn = 1;
if ~exist('stimlen','var') || isempty(stimlen)
    stimlen = str2double(regexp(spk.stimlength, '\d{1,3}(?=min)','match','once'));
end
dft = 1;
dff = 5;
stimfolder = 'I:\Ripple_Noise\downsampled_for_MID';

stimstr = ne_get_ripple_noise_stimulus(stimfolder, rn, dft, dff, stimlen);

[nespktrain.evoked_spktrain, nespktrain.evoked_edges] = ne_create_spktrain_from_stim_mat(spk, stimstr.stimulus, trigger);

sponend = (trigger(1) - 1)/20000*1000;
spon_edges = 0:0.5:sponend;
spon_spktrain = zeros(length(spk.spk), length(spon_edges)-1);
spiketimes = {spk.spk.spiketimes};

for j = 1:length(spk.spk)
    [spon_spktrain(j,:), ~] = histcounts(spiketimes{j}, spon_edges);
end

nespktrain.spon_edges = spon_edges;
nespktrain.spon_spktrain = spon_spktrain;

end

