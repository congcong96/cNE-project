function [nedata] = ne_calc_NE_for_spon_evoked_activity(spktrain, DF, edges, boundaryidx, position)
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

spontrain = spktrain(:, 1:boundaryidx - 1);
evokedtrain = spktrain(:, boundaryidx:end);

sponspkmat = downsample_spiketrain(spontrain, DF);
evokedspkmat = downsample_spiketrain(evokedtrain, DF);

spkmat = [sponspkmat evokedspkmat];

clear('spontrain','evokedtrain','sponedges','evokededges')

sponedges = edges(1:boundaryidx);
evokededges = edges(boundaryidx+1:end);

edges = [sponedges(1:DF:end) evokededges(1:DF:end)];

maxdiff = max(unique(diff(edges)));
edges = [edges edges(end)+maxdiff];

boundaryidx = size(sponspkmat,2) + 1;


% Find cell assemblies
fprintf('\nDetecting cNEs for whole spiketrain\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spkmat);

fprintf('\nDetecting cNEs for spontaneous spiketrain\n');
[sponPatterns, sponActivities] = ca_detect_cell_assemblies_data(sponspkmat);

fprintf('\nDetecting cNEs for evoked spiketrain\n');
[evokedPatterns, evokedActivities] = ca_detect_cell_assemblies_data(evokedspkmat);




% Assign data to struct for output argument
nedata.spktrain = spkmat;
nedata.df = DF;
nedata.position = position;
nedata.total_patterns = Patterns;
nedata.total_activities = Activities;
nedata.spon_patterns = sponPatterns;
nedata.spon_activities = sponActivities;
nedata.evoked_patterns = evokedPatterns;
nedata.evoked_activities = evokedActivities;
nedata.edges = edges;
nedata.boundaryidx = boundaryidx;


return;




% Code for significance testing, if we want it.
% nreps = 4;
% pval = 0.05;
% for i = 1:size(stamat,1)
%     [sta_sig, siglevel, rand_dist] = ...
%         ca_sig_sta_from_stim_obs_resp(stamat(i,:), resp(i,:), stim, nreps, pval);
%     figure;
%     subplot(1,2,1);
%     imagesc(reshape(stamat(i,:),nf,nlags));
%     subplot(1,2,2);
%     imagesc(reshape(sta_sig,nf,nlags));
% end % (for i)




















