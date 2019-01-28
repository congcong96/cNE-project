function [nedata] = ne_calc_rn_rnrep_neuronal_ensemble_sta(longspktrain, repspktrain, stim_mat, FsDVD, DF, position)
% ca_calc_cell_assembly_sta Estimate STAs from PCA/ICA cell assembly analysis
% 
%     [cadata] = ca_calc_cell_assembly_sta(spk, stimulus, trigger, FsDVD, totalDF)
%     --------------------------------------------------------------------
%     spk : struct array of spike times. Obtaines from saveSpikeSort.m
%
%     stimulus : ripple stimulus in matrix format
%
%     trigger : trigger for ripple stimulus
%
%     FsDVD : sampling of DVD Audio system
%
%     DF : downsampling factor (1 is 0.5 ms, 2 is 1 ms, etc.)
%
%     cadata : struct holding calculations. Has the form:
%         cadata.spktrain = binned spike train, one row per neuron
%         cadata.fsdvd = sampling rate of sound stimulus
%         cadata.df = total downsampling factor of ripple noise envelope
%         cadata.position = recording depth. position(i) -> spktrain(i,:)
%         cadata.Patterns = cell assemblies
%         cadata.Activities = time course of cell assembly activity; one row per assembly
%         cadata.nf = number of frequencies in STA
%         cadata.nlags = number of time lags in STA
%         cadata.stamat = matrix of spike train STAs.
%         cadata.ca_stamat = matrix of cell assembly STAs.
%
%   updated 7/20/16 JS


% Make sure there are 6 input arguments.
% narginchk(6,6);

longspkmat = downsample_spiketrain(longspktrain, DF);
repspkmat = downsample_spiketrain(repspktrain, DF);
stimulus = stim_mat(:,1:DF:end);
boundaryidx = size(longspkmat, 2) + 1;

spkmat = [longspkmat repspkmat];

% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spkmat);
ca_plot_cell_assemblies_data(spkmat, Patterns, Activities, position);


% Get stimulus-response matrices for later STA calculations
if DF <= 5 
    nlags = 40;
else
    nlags = 20;
end

if DF > 10
    sta_spktrain = downsample_spiketrain(longspktrain, 10);
    sta_stimulus = stim_mat(:,1:10:end);
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(sta_stimulus, sta_spktrain, nlags);
else
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(stimulus, spkmat, nlags);
end


% Calculate STAs in one block for each neuron.

stamat = resp * stim; % one sta per row
[nf, ~] = size(stimulus);
nlags = size(stamat,2) / nf; % Get # time lags from sta and stimulus

% Assign data to struct for output argument
nedata.spktrain = spkmat;
nedata.fsdvd = FsDVD;
nedata.df = DF;
nedata.position = position;
nedata.Patterns = Patterns;
nedata.Activities = Activities;
nedata.nf = nf;
nedata.nlags = nlags;
nedata.stamat = stamat;
nedata.boundaryidx = boundaryidx;


if DF > 10
    nedata.sta_spktrain = sta_spktrain;
end



return;


















