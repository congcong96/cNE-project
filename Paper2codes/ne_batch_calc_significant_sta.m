function ne_batch_calc_significant_sta(nefiles)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

pval = 95;
nreps = 100;

for i = 1:length(nefiles)
    
    clc;
    fprintf('\nProcessing significant STAs for %s (file %d of %d)...\n', nefiles{i}, i, length(nefiles))

    load(nefiles{i}, 'exp_site_nedata')
    
    nedata = exp_site_nedata.nedata;
    nlags = nedata.nlags;
    neuron_stamat = nedata.stamat;    
    NE_stamat = nedata.NE_stamat;  
    all_stamat = [NE_stamat; neuron_stamat];
    
    neuron_locator = nedata.sta_spktrain;
    NE_locator = nedata.sta_NEtrain;
    all_locator = [NE_locator; neuron_locator];
    
    dft = 10;
    dff = 5;
    rn = exp_site_nedata.stim; % stim as a string
    rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
    rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
    stimlength = exp_site_nedata.stimlength;

    [stimstr] = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlength);
    
    [sig_stamat, ptd, moransI, rel_idx, siglevel] = ...
        ne_calc_sig_sta_statistics(all_stamat, all_locator,...
        stimstr.stimulus, nreps, nlags, stimlength, pval);

    exp_site_nedata.nedata.sig_NE_stamat = sig_stamat(1:size(NE_stamat,1), :);
    exp_site_nedata.nedata.sig_stamat = sig_stamat(size(NE_stamat,1)+1:end, :);
    
    exp_site_nedata.nedata.sta_NE_siglevel = siglevel(1:size(NE_stamat,1), :);
    exp_site_nedata.nedata.sta_siglevel = siglevel(size(NE_stamat,1)+1:end, :);
    
    % peak-trough difference
    shuffled_stats.NE_ptd = ptd(1:size(NE_stamat,1), :); 
    shuffled_stats.neuron_ptd = ptd(size(NE_stamat,1)+1:end, :);
    
    % Moran's I
    shuffled_stats.NE_moransI = moransI(1:size(NE_stamat,1), :);
    shuffled_stats.neuron_moransI = moransI(size(NE_stamat,1)+1:end, :);
    
    % STA reliability index
    shuffled_stats.NE_reliability_idx = rel_idx(1:size(NE_stamat,1),:);
    shuffled_stats.neuron_reliability_idx = rel_idx(size(NE_stamat,1)+1:end, :);
    
    exp_site_nedata.nedata.shuffled_stats = shuffled_stats;
    
    save(nefiles{i}, 'exp_site_nedata', '-append');
    clear('exp_site_nedata')

end

