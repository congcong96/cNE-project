function NEneuroninfo = ne_calc_ne_neuron_sta_info(exp_site_nedata, sigcheckopt, NEneuroninput)

% NEneuroninput: For updating of existing NEneuroninfo struct arrays

if ~exist('NEneuroninput','var')
    NEneuroninput = [];
end

fraction = [90 95 97.5 99 100];
spike_count_threshold = 300; % only calc info if each subset has this many spikes
nsamples = 50; % number of times to subsample to make num of spks equal in all conditions
niter = 20; % number of times to calc info of frac of sample

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;
% neact = nedata.Activities;
% nethresh = nedata.NEthresh;
members = nedata.NEmembers;

if all(cellfun('isempty', members))
    NEneuroninfo = [];
    return
end

% if size(nethresh,1) > 1
%     alpha = nedata.NEthresh_alpha;
%     idx = NEthreshalpha == alpha;
%     if ~any(idx)
%         error('Choose valid NEthresh_alpha')
%     end
%     nethresh = nethresh(idx,:);
% end


if isfield(exp_site_nedata, 'stimlength')
    stimlen = exp_site_nedata.stimlength;
else
    stimlen = 10;
end


if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
    spktrain = nedata.spktrain;
else
    dft = 10;
    spktrain = nedata.sta_spktrain;
end

% get neurons part of at least one cNE
% neunespks = cellfun(@(x) logical(spktrain(x,:)), members, 'UniformOutput', 0);
% spkcount = cellfun(@(x) sum(x, 2), neunespks, 'UniformOutput', 0);

% get cNE events
NEraster = nedata.sta_NEtrain;

if sigcheckopt
    
    NEnum = find(nedata.sig_NE_sta);
    signeuron = find(nedata.sig_neuron_sta);

    if ~any(NEnum)
        NEneuroninfo = [];
        return
    end
    
    % get significant cNEs and their event trains only
    NEraster = NEraster(NEnum, :);
    members = members(NEnum);

    % get significant neurons and their spike trains within significant cNEs only
    members = cellfun(@(x) intersect(x, signeuron), members, 'UniformOutput', 0); 
    neunespks = cellfun(@(x) logical(spktrain(x,:)), members, 'UniformOutput', 0);
    
else
    NEnum = 1:length(members);
    neunespks = cellfun(@(x) logical(spktrain(x,:)), members, 'UniformOutput', 0);
end

NEeventcount = sum(NEraster,2);

% initialize info struct array
totalneurons = sum(cellfun(@length, members));

if totalneurons == 0
    NEneuroninfo = [];
    return
end

spkcount = cellfun(@(x) sum(x, 2), neunespks, 'UniformOutput', 0);

if isempty(NEneuroninput)
    NEneuroninfo(totalneurons).neuron = [];
else
    NEneuroninfo(totalneurons) = NEneuroninput(1);
    allneurons = [NEneuroninput.neuron];
    allNEs = [NEneuroninput.NE];
end


% get stimulus
dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';

stimstr = ca_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlen);


c = 1;

for j = 1:length(spkcount)
    
    fprintf('NE #%d vs neuron info: \n', NEnum(j));
    
    for k = 1:length(spkcount{j})
        
        if ~isempty(NEneuroninput)
            index = find(members{j}(k) == allneurons & NEnum(j) == allNEs, 1);
            if ~isempty(index)
                NEneuroninfo(c) = NEneuroninput(index);
                fprintf('\n');
                fprintf('Information for NE-neuron pairs:\n');
                fprintf('NE %d: %.3f bits/spk\n', NEnum(j), NEneuroninfo(c).NE_info_extrap);
                fprintf('Neuron %d: %.3f bits/spk\n', members{j}(k), NEneuroninfo(c).neuron_info_extrap);
                fprintf('\n');
                c = c+1;
                continue
            end
        end        
                               
        NEneuroninfo(c).neuron = members{j}(k);
        NEneuroninfo(c).NE = NEnum(j);
        
        numcheck = min([NEeventcount(j) spkcount{j}(k)]);

        if numcheck > spike_count_threshold

            %get minimum number of spikes (95% of the smallest group)
            min_spikes = round(0.95*numcheck);

            %initialize cell array to store all info values                
            NE_ifraction = cell(nsamples,1);
            neuron_ifraction = cell(nsamples,1);  

            for i = 1:nsamples

                samp_NE = sub_sample_spktrain(NEraster(j,:), NEeventcount(j) - min_spikes);
                samp_neuron = sub_sample_spktrain(neunespks{j}(k,:), spkcount{j}(k) - min_spikes);                

                sta_NE = calc_single_sta_from_locator_stimulus(samp_NE, stimstr.stimulus, nlags);
                sta_neuron = calc_single_sta_from_locator_stimulus(samp_neuron, stimstr.stimulus, nlags);

                [xprior_NE, xposterior_NE] = ne_sta_stimulus_projection(sta_NE, samp_NE, stimstr.stimulus);
                [xprior_neuron, xposterior_neuron] = ne_sta_stimulus_projection(sta_neuron, samp_neuron, stimstr.stimulus);                

                NE_ifraction{i} = ca_subset_info_from_data_fraction2(xprior_NE, xposterior_NE, fraction, niter);
                neuron_ifraction{i} = ca_subset_info_from_data_fraction2(xprior_neuron, xposterior_neuron, fraction, niter);

            end

            % Get mean/std of information for the data fractions

            [NE_info_frac_mn, NE_info_frac_std] = ...
                info_from_data_fraction_mean_std(fraction, NE_ifraction);

            [neuron_info_frac_mn, neuron_info_frac_std] = ...
                info_from_data_fraction_mean_std(fraction, neuron_ifraction);


            % Extrapolate the information values to get the final value for each spike train type
            [NE_info_extrap] = info_extrapolate_from_mn_std(fraction, NE_info_frac_mn);
            [neuron_info_extrap] = info_extrapolate_from_mn_std(fraction, neuron_info_frac_mn);


            fprintf('\n');
            fprintf('Information for NE-neuron pairs:\n');
            fprintf('NE %d: %.3f bits/spk\n', NEnum(j), NE_info_extrap);
            fprintf('Neuron %d: %.3f bits/spk\n', members{j}(k), neuron_info_extrap);
            fprintf('\n');


            %save data

            NEneuroninfo(c).fraction = fraction;
            NEneuroninfo(c).min_spikes = min_spikes;
            NEneuroninfo(c).NEmembers = nedata.NEmembers{NEnum(j)};

            NEneuroninfo(c).NE_info_frac_mn = NE_info_frac_mn;
            NEneuroninfo(c).NE_info_frac_std = NE_info_frac_std;
            NEneuroninfo(c).NE_info_extrap = NE_info_extrap;

            NEneuroninfo(c).neuron_info_frac_mn = neuron_info_frac_mn;
            NEneuroninfo(c).neuron_info_frac_std = neuron_info_frac_std;
            NEneuroninfo(c).neuron_info_extrap = neuron_info_extrap;

            c = c+1;

        else
            %if one or more groups did not cross threshold, save empty arrays
            NEneuroninfo(c).fraction = [];
            NEneuroninfo(c).min_spikes = [];
            NEneuroninfo(c).NEmembers = [];

            NEneuroninfo(c).NE_info_frac_mn = [];
            NEneuroninfo(c).NE_info_frac_std = [];
            NEneuroninfo(c).NE_info_extrap = [];

            NEneuroninfo(c).neuron_info_frac_mn = [];
            NEneuroninfo(c).neuron_info_frac_std = [];
            NEneuroninfo(c).neuron_info_extrap = [];

            c = c+1;
        end

    end       
end
       
end
                
                


