function ne_batch_split_spon_evoked_spkstruct(spkfiles, stimlen)

if ~exist('stimlen','var')
    stimlen = [];
end

for i = 1:length(spkfiles)
    load(spkfiles{i}, 'spk','trigger')
    nespktrain = ne_split_spon_evoked_spkstruct(spk, trigger, stimlen);
    save(spkfiles{i}, 'spk', 'trigger')
    save(spkfiles{i}, '-struct', 'nespktrain', '-append'); 
end