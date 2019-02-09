function neuronNEinfo = ne_calc_NEsta_vs_neuronsta_info(exp_site_nedata, NEopt)

% we can do this 2 ways: project IC weights back onto the sta_spiketrain,
% or interpolate NE activity. 

% initialize standard values
spkcountthresh = 300;
nedata = exp_site_nedata.nedata;
stadf = 10;
nlags = 50;
fraction = [90 95 97.5 99 100];
nsamples = 10;

switch NEopt
    case 'interpolate' % interpolate NEactivities
        activities = nedata.Activities;
        NEact = zeros(size(activities,1), size(activities,2) * 2);
        for k = 1:size(activities,1)
            NEact(k,:) = interp(activities(k,:), 2);
        end
    case 'project' % project IC weights back onto sta_spiketrain
        NEact = assembly_activity(nedata.Patterns, nedata.sta_spktrain);
    otherwise
        error('NEopt must be ''interpolate'' or ''project''!')
end
        
spktrain = nedata.sta_spktrain;

% make spktrain and NEact similar in length
if size(spktrain,2) ~= size(NEact,2)
    diff = size(spktrain,2) - size(NEact,2);
    if diff > 0 
        spktrain = spktrain(:,1:end-diff);
    else
        NEact = NEact(:,1:end-diff);
    end
end

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

% make stimulus matrix the same length as everything else
if size(stimulus,2) ~= size(NEact,2)
    diff = size(stimulus,2) - size(NEact,2);
    stimulus = stimulus(:,1:end-diff);
end

% get NE members
NEmembers = nedata.NEmembers;

%initialize neuronNEinfo
total = sum(cellfun('length',NEmembers));
neuronNEinfo(total).exp = [];

c = 1;

for i = 1:length(NEmembers)
    
    NEtrain = NEact(i,:); 
    
    for j = 1:length(NEmembers{i})
        
        neuronNEinfo(c).exp = exp_site_nedata.exp;
        neuronNEinfo(c).site = exp_site_nedata.site;
        neuronNEinfo(c).stim = exp_site_nedata.stim;
        
        % get neuron rasters and min_events
        neuronraster = logical(spktrain(NEmembers{i}(j),:));
        neuronsum = sum(neuronraster);        
        min_events = round(0.95 * neuronsum);
        
        if min_events > spkcountthresh
            [NEraster, ~] = ne_calc_fixed_number_NE_rasters(NEtrain, neuronsum);
            
            for ii = 1:nsamples
                
                % Sub sample neuron and NE rasters
                samp_neuron = sub_sample_spktrain(neuronraster, neuronsum - min_events);
                samp_NE = sub_sample_spktrain(NEraster, neuronsum - min_events);

                % Estimate STAs for each group
                sta_neuron = ca_sta_from_locator_stimulus(samp_neuron, stimulus, nlags);
                sta_NE = ca_sta_from_locator_stimulus(samp_NE, stimulus, nlags);
                
                % Get all the projection values; prior means without regard to a
                % spike; posterior means for a spike
                [xprior_neuron, xposterior_neuron] = ne_sta_stimulus_projection(sta_neuron, samp_neuron, stimulus);
                [xprior_NE, xposterior_NE] = ne_sta_stimulus_projection(sta_NE, samp_NE, stimulus);


                % Calculate info for different data fractions             
                neuron_ifraction{ii} = ca_subset_info_from_data_fraction2(xprior_neuron, xposterior_neuron, fraction);
                NE_ifraction{ii} = ca_subset_info_from_data_fraction2(xprior_NE, xposterior_NE, fraction);

            end



            % Get mean/std of information for the data fractions

            [neuron_info_frac_mn, neuron_info_frac_std] = ...
                info_from_data_fraction_mean_std(fraction, neuron_ifraction);

            [NE_info_frac_mn, NE_info_frac_std] = ...
                info_from_data_fraction_mean_std(fraction, NE_ifraction);


            % Extrapolate the information values to get the final value for each spike train type
            [neuron_info_extrap] = info_extrapolate_from_mn_std(fraction, neuron_info_frac_mn);
            [NE_info_extrap] = info_extrapolate_from_mn_std(fraction, NE_info_frac_mn);


            fprintf('\n');
            fprintf('Information for neuron-NE pairs:\n');
            fprintf('neuron #%d: %.3f bits/spk\n', NEmembers{i}(j), neuron_info_extrap);
            fprintf('NE #%d: %.3f bits/spk\n', i, NE_info_extrap);
            
            
            % Data from raw STA
            neuronNEinfo(c).fraction = fraction;
            neuronNEinfo(c).NE = i;
            neuronNEinfo(c).neuron = NEmembers{i}(j);
            neuronNEinfo(c).min_spikes = min_events;

            neuronNEinfo(c).neuron_info_frac_mn = neuron_info_frac_mn;
            neuronNEinfo(c).neuron_info_frac_std = neuron_info_frac_std;
            neuronNEinfo(c).neuron_info_extrap = neuron_info_extrap;

            neuronNEinfo(c).NE_info_frac_mn = NE_info_frac_mn;
            neuronNEinfo(c).NE_info_frac_std = NE_info_frac_std;
            neuronNEinfo(c).NE_info_extrap = NE_info_extrap;
            
        else
            neuronNEinfo(c).fraction = fraction;
            neuronNEinfo(c).NE = i;
            neuronNEinfo(c).neuron = NEmembers{i}(j);
            neuronNEinfo(c).min_spikes = min_events;

            neuronNEinfo(c).neuron_info_frac_mn = [];
            neuronNEinfo(c).neuron_info_frac_std = [];
            neuronNEinfo(c).neuron_info_extrap = [];

            neuronNEinfo(c).NE_info_frac_mn = [];
            neuronNEinfo(c).NE_info_frac_std = [];
            neuronNEinfo(c).NE_info_extrap = [];
            
            
        end
        c = c+1;
    end
end

                
                     

        
        
        