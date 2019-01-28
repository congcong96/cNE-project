function trialspk = ne_calc_rnrep_trialspk(spktimes, trigger, offset)

% Separates spiketimes from rnrep data into their respective trials.
% 
%   Inputs:
%       spktimes: time of spikes in ms.
%       trigger: time of trigger in sample number.
%       offset: 'exclusion zone' to eliminate onset activity in ms.
% 
%   Outputs:
%       trialspk: cell array of spiketimes of each trial (in each cell)
% 
%   Jermyn See, updated 7/17/15.


stimlength = 4700;  %ms, actually 4669, but rounded up to include post-stim responses
% offset = 200;      %ms

trigms = trigger/20000*1000;

trialstart = trigms + offset;
trialend = trigms + stimlength;

trialspk = arrayfun(@(x,y,z) spktimes(spktimes >= x & spktimes <= y) - z,...
    trialstart, trialend, trigms, 'UniformOutput',0);
    


