function [FSI_ne, FSI_neuron] = ne_calc_NE_neuron_FSI(exp_site_nedata)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;
NEmembers = nedata.NEmembers;
NEact = nedata.Activities;
spktrain = nedata.sta_spktrain;
NEthresh = nedata.NEthresh;


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
FSI_ne = zeros(size(sta,1),1);

binedges = -6:0.1:6;
ideal_cdf = zeros(1, length(binedges)-1);
ideal_cdf(end) = 1;

for i = 1:size(sta,1)

    [xprior, xposterior] = ne_sta_stimulus_projection(reshape(sta(i,:),nf, nlags), NEtrain(i,:), stimstr.stimulus);
    xprior(1:nlags-1) = []; % remove zeros at the start of the stimulus
    
    xmn = mean(xprior);
    xstd = std(xprior);
    xprior_norm = (xprior - xmn) ./ xstd;
    xposterior_norm = (xposterior - xmn) ./ xstd;
    
    xprior_cdf = cumsum(histcounts(xprior_norm, binedges, 'Normalization','probability'));
    xposterior_cdf = cumsum(histcounts(xposterior_norm, binedges, 'Normalization','probability'));
    
    
    FSI_ne(i) = sum(xprior_cdf-xposterior_cdf)/sum(xprior_cdf-ideal_cdf);
   


%     [xposterior_cdf, cdfbins] = ecdf(xposterior_norm);
% 
%     cdfidx = arrayfun(@(x) find(cdfbins<=x, 1,'last'), FSIthresh,'UniformOutput',0);
%     cdfidx(cellfun('isempty',cdfidx)) = {1};
%     cdfidx = cell2mat(cdfidx);
%     FSI_ne(i,:) = xposterior_cdf(cdfidx);

    
%     FSI_ne(i) = sum(xposterior_norm >= FSIthresh) / length(xposterior_norm);
%     NFSI_ne(i) = sum(xposterior_norm <= NFSIthresh) /length(xposterior_norm);
%     xposterior_shuff_norm{i} = (xposteriorshuff(:) - xmn) ./ xstd;
    
%     figure; hold on
%     histogram(xprior_norm, xbin_edges, 'Normalization', 'probability')
%     histogram(xposterior_norm, xbin_edges, 'Normalization','probability')
%     histogram(xposterior_shuff_norm{i}, xbin_edges{i}, 'Normalization','probability')
%     legend('PV shuffled', 'PV data')

end

neuronsta = nedata.stamat;
sta_spktrain = nedata.sta_spktrain;
FSI_neuron = zeros(size(neuronsta,1), 1);


for i = 1:size(neuronsta,1)
    [xprior, xposterior] = ne_sta_stimulus_projection(reshape(neuronsta(i,:),nf, nlags), sta_spktrain(i,:), stimstr.stimulus);
    xprior(1:nlags-1) = []; % remove zeros at the start of the stimulus
    
        
    xmn = mean(xprior);
    xstd = std(xprior);
    xprior_norm = (xprior - xmn) ./ xstd;
    xposterior_norm = (xposterior - xmn) ./ xstd;

    xprior_cdf = cumsum(histcounts(xprior_norm, binedges, 'Normalization','probability'));
    xposterior_cdf = cumsum(histcounts(xposterior_norm, binedges, 'Normalization','probability'));
    
    FSI_neuron(i) = sum(xprior_cdf-xposterior_cdf)/sum(xprior_cdf-ideal_cdf);
   
%     xbin_edges = linspace(min(xprior_norm{i}), max(xprior_norm{i}), numbins+1);

%     xposterior_norm = (xposterior - xmn) ./ xstd;
%     [xposterior_cdf, cdfbins] = ecdf(xposterior_norm);
% 
%     cdfidx = arrayfun(@(x) find(cdfbins<=x, 1,'last'), FSIthresh,'UniformOutput',0);
%     cdfidx(cellfun('isempty',cdfidx)) = {1};
%     cdfidx = cell2mat(cdfidx);
%     FSI_neuron(i,:) = xposterior_cdf(cdfidx);
%     
end
% 
% numNEmembers = sum(cellfun('length', NEmembers));
% FSI_paired = zeros(numNEmembers, 2);
% c = 1;
% 
% for i = 1:length(NEmembers)
%     
%     for j = 1:length(NEmembers{i})
%         
%         FSI_paired(c,1) = FSI_ne(i);
%         FSI_paired(c,2) = FSI_neuron(NEmembers{i}(j));
%         c = c+1;
%     end
% end