function pseudoNEinfo = ne_calc_NEsta_vs_pseudoNEsta_info(exp_site_nedata)

% we can do this 2 ways: project IC weights back onto the sta_spiketrain,
% or interpolate NE activity. 

% initialize standard values
spkcountthresh = 300;
nedata = exp_site_nedata.nedata;
stadf = 10;
nlags = 50;
fraction = [90 95 97.5 99 100];
nsamples = 10;
% spktrain = nedata.spktrain;
spktrain = nedata.sta_spktrain;
randiter = 50;
NEmembers = nedata.NEmembers;
NEsize = cellfun('length', NEmembers);
% NEthresh = nedata.NEthresh;

% % interpolate NEactivities for STA calculations
% activities = nedata.Activities;
% NEact = zeros(size(activities,1), size(activities,2) * 2);
% NEraster = zeros(size(activities,1), size(activities,2) * 2);
% 
% for k = 1:size(activities,1)
%     NEact(k,:) = interp(activities(k,:), 2);
%     NEraster(k,:) = NEact(k,:) >= NEthresh(k);
% end

NEraster = cell2mat(cellfun(@(x) sum(spktrain(x,:),1), NEmembers, 'UniformOutput' ,0)) > 0;

numevents = sum(NEraster,2);

% pNEact = cell(randiter, 1);
pNEmembers = cell(randiter, 1);
pNEraster = cell(randiter, 1);

% for jj = 1:randiter
%     pICweights = ne_shuffle_IC_weights(exp_site_nedata, 'random');
%     temp = assembly_activity(pICweights, spktrain);
%     pNEthresh = ne_calc_NE_act_thresholds(exp_site_nedata,'circular', 10, 99.9, pICweights);
%     
%     for j = 1:size(temp,1)
%         pNEact{jj}(j,:) = interp(temp(j,:), 2);
%         pNEraster{jj}(j,:) = pNEact{jj}(j,:) >= pNEthresh(j);
%     end 
%     
% 
% end

for j = 1:randiter
    pNEmembers{j} = arrayfun(@(x) randsample(1:size(spktrain,1), x), NEsize, 'UniformOutput', 0);
    pNEraster{j} = cell2mat(cellfun(@(x) sum(spktrain(x,:),1), pNEmembers{j}, 'UniformOutput' ,0)) > 0;
end

pnumevents = cell2mat(cellfun(@(x) sum(x,2), pNEraster', 'UniformOutput', 0))';
        
% get stimulus matrix
drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);

stimtype = exp_site_nedata.stim;
stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stimtype);
stimfile = fullfile(stimfolder,stimfilename);
stimfile = dir(stimfile);
load(stimfile.name);

stimulus = stim_mat(:,1:stadf:end);
clear('stim_mat');

% make stimulus matrix same size as everything else
if size(stimulus,2) ~= size(NEraster,2)
    diff = size(stimulus,2) - size(NEraster,2);
    stimulus = stimulus(:,1:end-diff);
end

pseudoNEinfo(size(NEraster,1)).exp = [];

for i = 1:size(NEraster, 1)
    
%     NEtrain = NEraster(i,:); 

    fprintf('\nCalculating NE #%d of %d...\n', i, size(NEraster,1))
       
    pseudoNEinfo(i).exp = exp_site_nedata.exp;
    pseudoNEinfo(i).site = exp_site_nedata.site;
    pseudoNEinfo(i).stim = exp_site_nedata.stim;

    % get min_events
    min_events = round(0.95 * min([numevents(i); pnumevents(:,i)]));

    if min_events > spkcountthresh
        
        fprintf('\nCalculating real NE info...')

        for ii = 1:nsamples

            % Sub sample NE rasters
            samp_NE = sub_sample_spktrain(NEraster(i,:), numevents(i) - min_events);
%             samp_NE = sub_sample_spktrain(NEraster(i,:), round(0.05*numevents(i)));

            % Estimate STAs
            sta_NE = ca_sta_from_locator_stimulus(samp_NE, stimulus, nlags);

            % Get all the projection values; prior means without regard to a
            % spike; posterior means for a spike
            [xprior_NE, xposterior_NE] = ca_sta_stimulus_projection(sta_NE, samp_NE, stimulus);


            % Calculate info for different data fractions             
            NE_ifraction{ii} = ca_subset_info_from_data_fraction2(xprior_NE, xposterior_NE, fraction);

        end



        % Get mean/std of information for the data fractions

        [NE_info_frac_mn, NE_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, NE_ifraction);


        % Extrapolate the information values to get the final value for each spike train type
        
        [NE_info_extrap] = info_extrapolate_from_mn_std(fraction, NE_info_frac_mn);
        
        pNE_info_frac_mn = cell(randiter,1);
        pNE_info_frac_std = cell(randiter,1);
        pNE_info_extrap_all = zeros(randiter,1);
        
        for kk = 1:randiter
            
            fprintf('\nCalculating pseudo NE info, iteration %d of %d...', kk, randiter)
            
            for iii = 1:nsamples

            % Sub sample NE rasters
            samp_pNE = sub_sample_spktrain(pNEraster{kk}(i,:), pnumevents(kk,i) - min_events);
%             samp_pNE = sub_sample_spktrain(pNEraster{kk}(i,:), round(0.05*pnumevents(kk,i)));

            % Estimate STAs
            sta_pNE = ca_sta_from_locator_stimulus(samp_pNE, stimulus, nlags);

            % Get all the projection values; prior means without regard to a
            % spike; posterior means for a spike
            [xprior_pNE, xposterior_pNE] = ca_sta_stimulus_projection(sta_pNE, samp_pNE, stimulus);


            % Calculate info for different data fractions             
            pNE_ifraction{kk}{iii} = ca_subset_info_from_data_fraction2(xprior_pNE, xposterior_pNE, fraction);

            end
            
        % Get mean/std of information for the data fractions

        [pNE_info_frac_mn{kk}, pNE_info_frac_std{kk}] = ...
            info_from_data_fraction_mean_std(fraction, pNE_ifraction{kk});


        % Extrapolate the information values to get the final value for each spike train type
        
        [pNE_info_extrap_all(kk)] = info_extrapolate_from_mn_std(fraction, pNE_info_frac_mn{kk});
            
        end
        
        pNE_info_extrap = mean(pNE_info_extrap_all);
        
        fprintf('\n\n');
        fprintf('Information for NE-pseudoNE pairs:\n');
        fprintf('NE #%d: %.3f bits/spk\n', i, NE_info_extrap);
        fprintf('pseudo-NE #%d: %.3f bits/spk\n', i, pNE_info_extrap);

        % Data from raw STA
        pseudoNEinfo(i).fraction = fraction;
        pseudoNEinfo(i).NE = i;
        pseudoNEinfo(i).NEevents = numevents(i);
        pseudoNEinfo(i).pNEevents = pnumevents(:,i);
        
        pseudoNEinfo(i).NE_info_frac_mn = NE_info_frac_mn;
        pseudoNEinfo(i).NE_info_frac_std = NE_info_frac_std;
        pseudoNEinfo(i).NE_info_extrap = NE_info_extrap;
        
        pseudoNEinfo(i).pNE_info_frac_mn = pNE_info_frac_mn;
        pseudoNEinfo(i).pNE_info_frac_std = pNE_info_frac_std;
        pseudoNEinfo(i).pNE_info_extrap_all = pNE_info_extrap_all;
        pseudoNEinfo(i).pNE_info_extrap = pNE_info_extrap;


    else
        pseudoNEinfo(i).fraction = [];
        pseudoNEinfo(i).NE = [];
        pseudoNEinfo(i).NEevents = [];
        pseudoNEinfo(i).pNEevents = [];
        
        pseudoNEinfo(i).NE_info_frac_mn = [];
        pseudoNEinfo(i).NE_info_frac_std = [];
        pseudoNEinfo(i).NE_info_extrap = [];
        
        pseudoNEinfo(i).pNE_info_frac_mn = [];
        pseudoNEinfo(i).pNE_info_frac_std = [];
        pseudoNEinfo(i).pNE_info_extrap_all = [];
        pseudoNEinfo(i).pNE_info_extrap = [];


    end


end

                
                     

        
        
        