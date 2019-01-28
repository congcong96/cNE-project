function NEneurons = ne_find_NE_pairs_or_groups(exp_site_nedata, num, numsample, coincheck)

% find pairs or groups of NE neurons
% 
% Inputs:
% 
% num: size of groups, e.g. enter 2 for pairs.
% 
% numsample: number of groups to calculate. If less than total, will
% subsample randomly.
% 
% coincheck: takes values of 0 or 1. Checks if groups fire together at
% least once.

if ~exist('numsample','var')
    numsample = [];
end
if ~exist('coincheck','var')
    coincheck = 0;
end

nedata = exp_site_nedata.nedata;
EM = nedata.NEmembers;

% only consider NEs that are larger than the number for comparison
nummem = cellfun(@length,EM);
idx = nummem >= num;

% get unique combinations
comb = cell2mat(cellfun(@(x) nchoosek(x, num), EM(idx),'UniformOutput', 0));
unicomb = unique(comb, 'rows');

if coincheck == 1
    coinratio = ne_calc_coincidence_within_spktrain(nedata.spktrain, unicomb);
    unicomb(coinratio == 0,:) = [];
end

%randomly samples the valid combinations for the required number, or if
%less than the required number, returns all and a warning that less than
%the required number was returned
if ~isempty(numsample) && numsample < size(unicomb,1)
    
    NEneurons = unicomb(randsample(1:size(unicomb,1),numsample),:);
    
elseif isempty(numsample)
    
    NEneurons = unicomb;
    
else %~isempty(numsample) && numsample > size(unicomb,1)
    
    NEneurons = unicomb;
    warning(['Not enough NE combinations to extract %d groups of %d, '...
        'returning %d groups instead'], numsample, num, size(NEneurons,1));
    
end
end

