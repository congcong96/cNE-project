function [corrprior, corrposterior] = ne_sta_stimulus_correlation(sta, locator, stimulus)
% ca_sta_stimulus_projection - projection for training and test data sets
% 
%     Calculates all projections onto filters over the entire stimulus
%     duration. Each filter was calculated from 1 of the 4 training sets.
% 
%     Input arguments:
% 
%     sta : receptive field matrix.
% 
%     locator : a vector that describes whether a spike occurred during each
%     of the stimulus trials. Values are >= 0.
% 
%     stimulus : the entire stimulus that was played to the neuron. It is
%     an nf x ntrials size matrix. nf represents the number for frequencies
%     in the stimulus. ntrials is the total number of time trials.
% 
%     Output arguments:
% 
%     xprior : all projections of the stimulus onto the sta
% 
%     xposterior : projections corresponding to a spike in locator
%
%     Updated 9/22/17 by JS to allow for multiple spikes per bin in the
%     xposterior.


assert(size(sta,1) == size(stimulus,1))

ntrials = size(stimulus,2);


% Project all trials onto filter to get the prior distribution.
% This distribution will include the training set projections and the
% test set projections. All trials will have a projection value
% associated with them.


% Find the prior projection values for all filters at the same time:
% corrprior = zeros(ntrials, 1);

[nf, nlags] = size(sta); % # frequencies, # time bins
stimvec = zeros(nf*nlags, ntrials - nlags + 1);


for i = nlags:ntrials
   temp = stimulus(:, i-nlags+1:i);
   stimvec(:,i) = temp(:);
end

corrprior = corr(stimvec, sta(:));

if ~isempty(locator)
    corrposterior = rude(locator, corrprior); % includes only posterior
    
%     if shuffleopt ~= 0
%         
%         shuffidx = randi(length(locator), shuffleopt, 1);
%         shuff_xposterior = zeros(shuffleopt, length(corrposterior));
%         
%         for j = 1:length(shuffidx)
%             shuffloc = circshift(locator, shuffidx(j));
%             shuff_xposterior(j,:) = rude(shuffloc, corrprior);
%         end
%         
%         varargout{1} = shuff_xposterior;
%     end
        

else
    corrposterior = [];
end

    
    
    


return;



