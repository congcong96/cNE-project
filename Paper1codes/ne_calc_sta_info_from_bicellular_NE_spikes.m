function [subsetinfo] = ne_calc_sta_info_from_bicellular_NE_spikes(exp_site_nedata, spktrain)
% ca_get_cell_assembly_sta_info  Info from cell assemblies
%
% [subsetinfo] = ca_calc_sta_info_from_subsets(exp_site_nedata)
% ------------------------------------------------------------------------
%
%

fraction = [90 95 97.5 99 100];
nsamples = 10;
spike_count_threshold = 100; % only calc info if each subset has this many spikes

drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);

DF = 10;
stim = exp_site_nedata.stim;
stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stim);
stimfile = fullfile(stimfolder,stimfilename);
stimfile = dir(stimfile);
load(stimfile.name);
stimulus = stim_mat(:,1:DF:end);

bicell = ne_get_bicellular_spktrain_subset(exp_site_nedata, spktrain);
spksubset = ne_match_bicell_NE_vs_nonNE(bicell);


nonNE_ifraction = cell(nsamples,1);
NE_ifraction = cell(nsamples,1);

subsetinfo(length(spksubset)).exp = [];

nlags = exp_site_nedata.nedata.nlags;

for i = 1:length(spksubset)

    subsetinfo(i).exp = exp_site_nedata.exp;
    subsetinfo(i).site = exp_site_nedata.site;
    subsetinfo(i).stim = exp_site_nedata.stim;
    subsetinfo(i).df = exp_site_nedata.df;
    subsetinfo(i).fs = exp_site_nedata.fs;
    
    subsetinfo(i).depth = exp_site_nedata.depth;
    subsetinfo(i).probetype = exp_site_nedata.probetype;
    subsetinfo(i).refneuron = spksubset(i).refneuron;
    subsetinfo(i).NEneuron = spksubset(i).NEneuron;
    subsetinfo(i).nonNEneuron = spksubset(i).nonNEneuron;
    subsetinfo(i).NEcomb = spksubset(i).NEcomb;
    subsetinfo(i).nonNEcomb = spksubset(i).nonNEcomb;
    


    % Get binned spike trains for each subset
    locator_wo_NE = spksubset(i).nonNEspktrain;
    locator_w_NE = spksubset(i).NEspktrain;

    fprintf('NE SUBSET INFO: %.0f of %.0f\n', i, length(spksubset));

    if ( sum(locator_wo_NE) > spike_count_threshold && ...
         sum(locator_w_NE) > spike_count_threshold )
     
       
        locator_wo_NE = locator_wo_NE > 0;
        locator_w_NE = locator_w_NE > 0;
        
        count_wo_NE = sum(locator_wo_NE);
        count_w_NE = sum(locator_w_NE);
        
        % Determine minimum number of samples/spikes. Datasets will be 
        % sampled based on this number.    
        min_spikes = min([count_wo_NE count_w_NE]);
        min_spikes = round(0.95*min_spikes);
        
        
        for j = 1:nsamples
            
            samp_wo_NE = sub_sample_spktrain(locator_wo_NE, count_wo_NE - min_spikes);
            samp_w_NE = sub_sample_spktrain(locator_w_NE, count_w_NE - min_spikes);

            % Estimate STAs for each subset
            sta_wo_NE = ca_sta_from_locator_stimulus(samp_wo_NE, stimulus, nlags);
            sta_w_NE = ca_sta_from_locator_stimulus(samp_w_NE, stimulus, nlags);

             % Plot STAs to make sure everything is working okay.
             % Uncomment below to check a neuron.
    %         figure;
    %         subplot(1,3,1); imagesc(sta_all);
    %         subplot(1,3,2); imagesc(sta_wo_NE);
    %         subplot(1,3,3); imagesc(sta_w_NE);
    %         pause;


            % Get all the projection values; prior means without regard to a
            % spike; posterior means for a spike
            [xprior_wo_NE, xposterior_wo_NE] = ca_sta_stimulus_projection(sta_wo_NE, samp_wo_NE, stimulus);
            [xprior_w_NE, xposterior_w_NE] = ca_sta_stimulus_projection(sta_w_NE, samp_w_NE, stimulus);


          % Calculate info for different data fractions
          % -----------------------------------------------------------
            nonNE_ifraction{j} = ca_subset_info_from_data_fraction2(xprior_wo_NE, xposterior_wo_NE, fraction);
            NE_ifraction{j} = ca_subset_info_from_data_fraction2(xprior_w_NE, xposterior_w_NE, fraction);
        
        end



            % Get mean/std of information for the data fractions

            [wo_NE_info_frac_mn, wo_NE_info_frac_std] = ...
                info_from_data_fraction_mean_std(fraction, nonNE_ifraction);

            [w_NE_info_frac_mn, w_NE_info_frac_std] = ...
                info_from_data_fraction_mean_std(fraction, NE_ifraction);


            % Extrapolate the information values to get the final value for each spike train type
            
            [wo_NE_info_extrap] = info_extrapolate_from_mn_std(fraction, wo_NE_info_frac_mn);
            [w_NE_info_extrap] = info_extrapolate_from_mn_std(fraction, w_NE_info_frac_mn);


            fprintf('\n');
            fprintf('Information for subsets:\n');
            fprintf('WO : %.3f bits/spk\n', wo_NE_info_extrap);
            fprintf('W  : %.3f bits/spk\n\n', w_NE_info_extrap);


            % Save the data
            %---------------------------------------------------------------

            % Data from raw STA
            subsetinfo(i).fraction = fraction; 
            subsetinfo(i).min_spikes = min_spikes;


            subsetinfo(i).w_NE_info_frac_mn = w_NE_info_frac_mn;
            subsetinfo(i).w_NE_info_frac_std = w_NE_info_frac_std;
            subsetinfo(i).w_NE_info_extrap = w_NE_info_extrap;

            subsetinfo(i).wo_NE_info_frac_mn = wo_NE_info_frac_mn;
            subsetinfo(i).wo_NE_info_frac_std = wo_NE_info_frac_std;
            subsetinfo(i).wo_NE_info_extrap = wo_NE_info_extrap;

        else % if there were no common spikes, then save as empty matrices

            % Save the data
            %---------------------------------------------------------------


            % Data from raw STA
            subsetinfo(i).fraction = [];
            subsetinfo(i).min_spikes = [];

       

            subsetinfo(i).w_NE_info_frac_mn = [];
            subsetinfo(i).w_NE_info_frac_std = [];
            subsetinfo(i).w_NE_info_extrap = [];


            subsetinfo(i).wo_NE_info_frac_mn = [];
            subsetinfo(i).wo_NE_info_frac_std = [];
            subsetinfo(i).wo_NE_info_extrap = [];


   end % if

%    if ( ~mod(i,25) )
%       save('temp-bicellular-sta-info-data-2', 'bcinfo', 'i');
%    elseif ( i == length(bcstr) )
%       save('temp-bicellular-sta-info-data-2', 'bcinfo', 'i');
%    end

end % (for i)

return;

















