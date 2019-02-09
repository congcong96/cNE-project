function info = ne_calc_NEsta_info_for_binsizes(fileid, stimtype, df, USopt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% df is empty vector if using all


if nargin == 2
    df = [];
end

% Initialize some standard values
stadf = 2;
nlags = 50;
fraction = [90 95 97.5 99 100];
nsamples = 10;

% Initialize stimulus
drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);

stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stimtype);
stimfile = fullfile(stimfolder,stimfilename);
stimfile = dir(stimfile);
load(stimfile.name);

stimulus = stim_mat(:,1:stadf:end);
clear('stim_mat');


switch USopt
    case 'interpolate'
        % Extract NE activity matrices for different binsizes
        NEact = ne_load_NE_act_or_pat_for_different_binsizes(fileid, stimtype, df, 'act');

        % Upsample NE activity matrices to common resolution (1 ms bins)
        actstruct = ne_interpolate_upsampled_activities_for_diff_binsizes(NEact);
    
    case 'project'
        % Extract NE pattern matrices for different binsizes
        ICweights = ne_load_NE_act_or_pat_for_different_binsizes(fileid, stimtype, df, 'pat');
        basefile = gfn(sprintf('%s*-%s-*cmb.mat', fileid, stimtype));
        load(basefile{1}, 'spktrain')
        spktrain = downsample_spiketrain(spktrain,stadf);
        
        actstruct = ne_project_upsampled_activities_for_diff_binsizes(ICweights, spktrain);
        
    case 'offset'
               
        actstruct = ne_batch_calc_upsampled_NEactivity(fileid, stimtype, df, stadf);
        fn = fieldnames(actstruct);
        
        for i = 1:length(fn)
            sizediff = size(actstruct.(fn{i}), 2) - size(stimulus,2);
            actstruct.(fn{i}) = actstruct.(fn{i})(:, 1:end-sizediff);
        end
        
end


% Threshold NE activity matrices to get logical
% [actbin, ~] = ne_threshold_NE_activity(actstruct);
[actbin, ~] = ne_threshold_NE_activity(actstruct, 'fixednum', 5000);

% Get minimum event count
eventcount = structfun(@(x) sum(x,2), actbin, 'UniformOutput', 0);
min_events = round(0.95*min(cell2mat(struct2cell(eventcount))));

fn = fieldnames(actbin);

for i = 1:length(fn)
    
    fprintf('\nProcessing %s neuronal ensembles...\n', fn{i})
    NEbin = actbin.(fn{i});
    
    % get stim matrix for specific bin size
    
    if size(actbin.(fn{i}),2) == size(stimulus, 2)
        tempstim = stimulus;
    else %if size(sta_act,2) > size(sta_spktrain, 2)
        difference = size(stimulus,2) - size(actbin.(fn{i}),2);
        tempstim = stimulus(:, 1:end-difference);
    end
    
    for j = 1:size(NEbin,1)        
        fprintf('Processing %d of %d neuronal ensembles...\n',j, size(NEbin,1))
        
        ifraction = cell(nsamples,1);
        
        for k = 1:nsamples
            fprintf('%d of %d samples\n', k, nsamples)
            
            % subsample eventrain
            samp = sub_sample_spktrain(NEbin(j,:), eventcount.(fn{i})(j) - min_events);
            
            % get sta of subsampled eventrain
            sta = ca_sta_from_locator_stimulus(samp, tempstim, nlags);
            
            % Get all the projection values; prior means without regard to a
            % spike; posterior means for a spike
            [xprior, xposterior] = ne_sta_stimulus_projection(sta, samp, tempstim);
            
            % Calculate info for different data fractions
            ifraction{k} = ca_subset_info_from_data_fraction2(xprior, xposterior, fraction);
            
        end
        
        % Get mean/std of information for the data fractions
        [info_frac_mn, info_frac_std] = info_from_data_fraction_mean_std(fraction, ifraction);
        
        % Extrapolate the information values to get the final value for each spike train type
        [info_extrap] = info_extrapolate_from_mn_std(fraction, info_frac_mn);
        
        % Print results
        fprintf('\n');
        fprintf('Information for NE #%d for %s:\n', j, fn{i});
        fprintf('All: %.3f bits/event\n', info_extrap);
        fprintf('\n');
        
        % Save results
        info.(fn{i})(j).fraction = fraction;
        info.(fn{i})(j).min_events = min_events;

        info.(fn{i})(j).info_frac_mn = info_frac_mn;
        info.(fn{i})(j).info_frac_std = info_frac_std;
        info.(fn{i})(j).info_extrap = info_extrap;

    end
    
end

        
       
        