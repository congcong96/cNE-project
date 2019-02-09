function FSI = ne_calc_NE_neuron_paired_FSI(exp_site_nedata)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;
NEmembers = nedata.NEmembers;
NEact = nedata.Activities;
spktrain = nedata.sta_spktrain;
NEthresh = nedata.NEthresh;
numNEmembers = sum(cellfun('length',NEmembers));
neuronsta = nedata.stamat;
sta_spktrain = nedata.sta_spktrain;

NEtrain = ne_upsample_NEact_using_member_neuron_activity(NEact, NEmembers, spktrain, NEthresh);

if ~isfield(exp_site_nedata, 'stimlength')
    stimlength = 10;
else
    stimlength = exp_site_nedata.stimlength;
end

if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
else
    dft = 10;
end

dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
nlags = nedata.nlags;
nf = nedata.nf;

[stimstr] = ca_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlength);

% FSIthresh = -3:0.5:3;

sta = quick_calc_sta(stimstr.stimulus, NEtrain, nlags, nf);


FSI = zeros(numNEmembers,1);

binedges = -6:0.1:6;
% ideal_cdf = zeros(1, length(binedges)-1);
% ideal_cdf(end) = 1;

c = 1;

for i = 1:size(sta,1)

    [xpriorNE, xposteriorNE] = ne_sta_stimulus_projection(reshape(sta(i,:),nf, nlags), NEtrain(i,:), stimstr.stimulus);
    xpriorNE(1:nlags-1) = []; % remove zeros at the start of the stimulus
    
    xmn = mean(xpriorNE);
    xstd = std(xpriorNE);
    xpriorNEnorm = (xpriorNE - xmn) ./ xstd;
    xposteriorNEnorm = (xposteriorNE - xmn) ./ xstd;
    
    xpriorNEcdf = cumsum(histcounts(xpriorNEnorm, binedges, 'Normalization','probability'));
    xposteriorNEcdf = cumsum(histcounts(xposteriorNEnorm, binedges, 'Normalization','probability'));
    
    for j = 1:length(NEmembers{i})
        
        [xpriorNeuron, xposteriorNeuron] = ne_sta_stimulus_projection...
            (reshape(neuronsta(NEmembers{i}(j),:),nf, nlags),...
            sta_spktrain(NEmembers{i}(j),:), stimstr.stimulus);
        
        xpriorNeuron(1:nlags-1) = []; % remove zeros at the start of the stimulus

        xmn = mean(xpriorNeuron);
        xstd = std(xpriorNeuron);
        xpriorNeuronnorm = (xpriorNeuron - xmn) ./ xstd;
        xposteriorNeuronnorm = (xposteriorNeuron - xmn) ./ xstd;

        xpriorNeuroncdf = cumsum(histcounts(xpriorNeuronnorm, binedges, 'Normalization','probability'));
        xposteriorNeuroncdf = cumsum(histcounts(xposteriorNeuronnorm, binedges, 'Normalization','probability'));    
    
        FSI(c) = sum(xpriorNEcdf - xposteriorNEcdf)/sum(xpriorNeuroncdf - xposteriorNeuroncdf);
        c = c+1;
    end
   
end

