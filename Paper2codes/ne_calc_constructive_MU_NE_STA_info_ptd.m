function con_sta_MI_PTD = ne_calc_constructive_MU_NE_STA_info_ptd(con_sta, NEfolder)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('NEfolder', 'var')
    NEfolder = 'I:\Cell_Assemblies\MS_NEs_paper2';
end

fraction = [90 95 97.5 99 100];
spike_count_threshold = 300; % only calc info if each subset has this many spikes
nsamples = 50; % number of times to subsample to make num of spks equal in all conditions
niter = 20; % number of times to calc info of frac of sample

con_sta_MI_PTD(length(con_sta)).filename = [];

for i = 1:length(con_sta)
    
    fprintf('\nConstructive cNE vs multi-unit STRF MI/PTD (comparison %d of %d):\n', i, length(con_sta))
    
    con_sta_MI_PTD(i).filename = con_sta(i).filename;
    con_sta_MI_PTD(i).NE = con_sta(i).NE;
    con_sta_MI_PTD(i).neurons = con_sta(i).neurons;
    
    load(fullfile(NEfolder, con_sta(i).filename), 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nlags = nedata.nlags;
    
    % get spike/event trains
    MUtrain = logical(sum(nedata.sta_spktrain(con_sta(i).neurons,:), 1));
    NEtrain = logical(nedata.sta_NEtrain(con_sta(i).NE,:));
    
    % get stimulus
    dft = exp_site_nedata.df;
    if dft > 10
        dft = 10;
    end    
    stimlen = exp_site_nedata.stimlength;
    dff = 5;
    rn = exp_site_nedata.stim; % stim as a string
    rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
    rnpath = 'I:\Ripple_Noise\downsampled_for_MID';

    stimstr = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlen);
    
    MUcount = sum(MUtrain);
    NEcount = sum(NEtrain);
    
    con_sta_MI_PTD(i).NE_count = NEcount;
    con_sta_MI_PTD(i).MU_count = MUcount;
    
    numcheck = min([MUcount NEcount]);    

    if numcheck > spike_count_threshold

        %get minimum number of spikes (95% of the smallest group)
        min_spikes = round(0.95*numcheck);

        %initialize cell array to store all info values                
        NE_ifraction = cell(nsamples,1);
        MU_ifraction = cell(nsamples,1);
        
        NE_ptd = zeros(1, nsamples);
        MU_ptd = zeros(1, nsamples);

        for j = 1:nsamples

            samp_NE = sub_sample_spktrain(NEtrain, NEcount - min_spikes);
            samp_MU = sub_sample_spktrain(MUtrain, MUcount - min_spikes);                

            sta_NE = calc_single_sta_from_locator_stimulus(samp_NE, stimstr.stimulus, nlags);
            sta_MU = calc_single_sta_from_locator_stimulus(samp_MU, stimstr.stimulus, nlags);
            
            NE_ptd(j) = max(sta_NE(:)) - min(sta_NE(:)) ./ min_spikes;
            MU_ptd(j) = max(sta_MU(:)) - min(sta_MU(:)) ./ min_spikes;

            [xprior_NE, xposterior_NE] = ne_sta_stimulus_projection(sta_NE, samp_NE, stimstr.stimulus);
            [xprior_MU, xposterior_MU] = ne_sta_stimulus_projection(sta_MU, samp_MU, stimstr.stimulus);                

            NE_ifraction{j} = ca_subset_info_from_data_fraction2(xprior_NE, xposterior_NE, fraction, niter);
            MU_ifraction{j} = ca_subset_info_from_data_fraction2(xprior_MU, xposterior_MU, fraction, niter);

        end
        
        fprintf('\n');
        fprintf('NE %d PTD: %.2f\n', con_sta(i).NE, mean(NE_ptd));
        fprintf('MU (neurons %s) PTD: %.2f\n', strjoin(cellstr(num2str(con_sta(i).neurons)),', '), mean(MU_ptd));

        % Get mean/std of information for the data fractions

        [NE_info_frac_mn, NE_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, NE_ifraction);

        [MU_info_frac_mn, MU_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, MU_ifraction);


        % Extrapolate the information values to get the final value for each spike train type
        [NE_info_extrap] = info_extrapolate_from_mn_std(fraction, NE_info_frac_mn);
        [MU_info_extrap] = info_extrapolate_from_mn_std(fraction, MU_info_frac_mn);


        fprintf('\n');
        fprintf('NE %d MI: %.3f bits/spk\n', con_sta(i).NE, NE_info_extrap);
        fprintf('MU (neurons %s) MI: %.3f bits/spk\n', strjoin(cellstr(num2str(con_sta(i).neurons)),', '), MU_info_extrap);
        fprintf('\n');


        %save data

        con_sta_MI_PTD(i).fraction = fraction;
        con_sta_MI_PTD(i).min_spikes = min_spikes;
        
        con_sta_MI_PTD(i).NE_ptd_vector = NE_ptd;
        con_sta_MI_PTD(i).NE_ptd = mean(NE_ptd);
        
        con_sta_MI_PTD(i).MU_ptd_vector = MU_ptd;
        con_sta_MI_PTD(i).MU_ptd = mean(MU_ptd);

        con_sta_MI_PTD(i).NE_info_frac_mn = NE_info_frac_mn;
        con_sta_MI_PTD(i).NE_info_frac_std = NE_info_frac_std;
        con_sta_MI_PTD(i).NE_info_extrap = NE_info_extrap;

        con_sta_MI_PTD(i).MU_info_frac_mn = MU_info_frac_mn;
        con_sta_MI_PTD(i).MU_info_frac_std = MU_info_frac_std;
        con_sta_MI_PTD(i).MU_info_extrap = MU_info_extrap;

    else
        %if one or more groups did not cross threshold, save empty arrays
        con_sta_MI_PTD(i).fraction = [];
        con_sta_MI_PTD(i).min_spikes = [];
        
        con_sta_MI_PTD(i).NE_ptd_vector = [];
        con_sta_MI_PTD(i).NE_ptd = [];
        
        con_sta_MI_PTD(i).MU_ptd_vector = [];
        con_sta_MI_PTD(i).MU_ptd = [];

        con_sta_MI_PTD(i).NE_info_frac_mn = [];
        con_sta_MI_PTD(i).NE_info_frac_std = [];
        con_sta_MI_PTD(i).NE_info_extrap = [];

        con_sta_MI_PTD(i).MU_info_frac_mn = [];
        con_sta_MI_PTD(i).MU_info_frac_std = [];
        con_sta_MI_PTD(i).MU_info_extrap = [];

    end

end      
  
    
end

