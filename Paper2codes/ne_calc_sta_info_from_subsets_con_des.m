function subsetinfo = ne_calc_sta_info_from_subsets_con_des(spksubset)

fraction = [90 95 97.5 99 100];
nsamples = 50;
niter = 20;
spike_count_threshold = 100; % only calc info if each subset has this many spikes

subsetinfo(length(spksubset)).exp = [];

for i = 1:length(spksubset)
    
    load(spksubset(i).filename, 'exp_site_nedata');
    nedata = exp_site_nedata.nedata;
    nlags = nedata.nlags;
    if exp_site_nedata.df <= 10
        dft = exp_site_nedata.df;
    else
        dft = 10;
    end
    dff = 5;
    rn = exp_site_nedata.stim; % stim as a string
    rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
    stimlen = exp_site_nedata.stimlength;
    rnpath = 'I:/Ripple_Noise/downsampled_for_MID';

    [stimstr] = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlen);
    stim_mat = stimstr.stimulus;
    

    subsetinfo(i).exp = exp_site_nedata.exp;
    subsetinfo(i).site = exp_site_nedata.site;
    subsetinfo(i).stim = exp_site_nedata.stim;
    subsetinfo(i).df = exp_site_nedata.df;
    subsetinfo(i).fs = exp_site_nedata.fs;

    subsetinfo(i).depth = exp_site_nedata.depth;
    subsetinfo(i).probetype = exp_site_nedata.probetype;
    subsetinfo(i).neuron = spksubset(i).neuron;
    subsetinfo(i).NE = spksubset(i).NE;
    
    subsetinfo(i).fraction = fraction;            

    % Get binned spike trains for each subset
    locator_all = spksubset(i).spktrain;
    locator_w_NE = spksubset(i).subset_spktrain;

    fprintf('\nNE SUBSET PTD/MI: %d of %d\n', i, length(spksubset));

    if ( sum(locator_all) > spike_count_threshold && ...
         sum(locator_w_NE) > spike_count_threshold )
     
        locator_all = locator_all > 0;
        locator_w_NE = locator_w_NE > 0;
        
        count_all = sum(locator_all);
        count_w_NE = sum(locator_w_NE);
        
        % Determine minimum number of samples/spikes. Datasets will be 
        % sampled based on this number.    
        min_spikes = min([count_all count_w_NE]);
        min_spikes = round(0.95*min_spikes);
        
        w_NE_ifraction = cell(nsamples,1);
        all_ifraction = cell(nsamples,1);        
        
        ptd_all = zeros(1, nsamples);
        ptd_w_NE = zeros(1, nsamples);
        
        for j = 1:nsamples
            
            samp_all = sub_sample_spktrain(locator_all, count_all - min_spikes);
            samp_w_NE = sub_sample_spktrain(locator_w_NE, count_w_NE - min_spikes);

            % Estimate STAs for each subset
            sta_all = calc_single_sta_from_locator_stimulus(samp_all, stim_mat, nlags);
            sta_w_NE = calc_single_sta_from_locator_stimulus(samp_w_NE, stim_mat, nlags);

            % Calculate STA PTD
            ptd_all(j) = max(sta_all(:)) - min(sta_all(:)) ./ min_spikes;
            ptd_w_NE(j) = max(sta_w_NE(:)) - min(sta_w_NE(:)) ./ min_spikes;


            % Get all the projection values; prior means without regard to a
            % spike; posterior means for a spike
            [xprior_all, xposterior_all] = ne_sta_stimulus_projection(sta_all, samp_all, stim_mat);
            [xprior_w_NE, xposterior_w_NE] = ne_sta_stimulus_projection(sta_w_NE, samp_w_NE, stim_mat);


          % Calculate info for different data fractions
          % -----------------------------------------------------------
            all_ifraction{j} = ca_subset_info_from_data_fraction2(xprior_all, xposterior_all, fraction, niter);
            w_NE_ifraction{j} = ca_subset_info_from_data_fraction2(xprior_w_NE, xposterior_w_NE, fraction, niter);
        
        end
        
        fprintf('\n');
        fprintf('PTD for subsets:\n');
        fprintf('All: %.2f\n', mean(ptd_all));
        fprintf('W  : %.2f\n', mean(ptd_w_NE));


        % Get mean/std of information for the data fractions

        [all_info_frac_mn, all_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, all_ifraction);

        [w_NE_info_frac_mn, w_NE_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, w_NE_ifraction);


        % Extrapolate the information values to get the final value for each spike train type
        [all_info_extrap] = info_extrapolate_from_mn_std(fraction, all_info_frac_mn);
        [w_NE_info_extrap] = info_extrapolate_from_mn_std(fraction, w_NE_info_frac_mn);


        fprintf('\n');
        fprintf('Information for subsets:\n');
        fprintf('All: %.3f bits/spk\n', all_info_extrap);
        fprintf('W  : %.3f bits/spk\n\n', w_NE_info_extrap);


        % Save the data
        %---------------------------------------------------------------

        % Data from raw STA
        subsetinfo(i).min_spikes = min_spikes;

        subsetinfo(i).all_ptd_vec = ptd_all;
        subsetinfo(i).all_ptd = mean(ptd_all);
        
        subsetinfo(i).w_NE_ptd_vec = ptd_w_NE;
        subsetinfo(i).w_NE_ptd = mean(ptd_w_NE);

        subsetinfo(i).all_info_frac_mn = all_info_frac_mn;
        subsetinfo(i).all_info_frac_std = all_info_frac_std;
        subsetinfo(i).all_info_extrap = all_info_extrap;

        subsetinfo(i).w_NE_info_frac_mn = w_NE_info_frac_mn;
        subsetinfo(i).w_NE_info_frac_std = w_NE_info_frac_std;
        subsetinfo(i).w_NE_info_extrap = w_NE_info_extrap;

    else % if there were no common spikes, then save as empty matrices

        % Save the data
        %---------------------------------------------------------------


        % Data from raw STA
        
        locator_all = locator_all > 0;
        locator_w_NE = locator_w_NE > 0;
        
        count_all = sum(locator_all);
        count_w_NE = sum(locator_w_NE);
        
        min_spikes = min([count_all count_w_NE]);
        min_spikes = round(0.95*min_spikes);
        
        subsetinfo(i).min_spikes = min_spikes;
        
        subsetinfo(i).all_ptd_vec = [];
        subsetinfo(i).all_ptd = [];
        
        subsetinfo(i).w_NE_ptd_vec = [];
        subsetinfo(i).w_NE_ptd = [];

        subsetinfo(i).all_info_frac_mn = [];
        subsetinfo(i).all_info_frac_std = [];
        subsetinfo(i).all_info_extrap = [];

        subsetinfo(i).w_NE_info_frac_mn = [];
        subsetinfo(i).w_NE_info_frac_std = [];
        subsetinfo(i).w_NE_info_extrap = [];


   end % if

%    if ( ~mod(i,25) )
%       save('temp-bicellular-sta-info-data-2', 'bcinfo', 'i');
%    elseif ( i == length(bcstr) )
%       save('temp-bicellular-sta-info-data-2', 'bcinfo', 'i');
%    end

end % (for i)

return;

















