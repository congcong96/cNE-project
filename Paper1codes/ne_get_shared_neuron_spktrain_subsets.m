function [spkstruct] = ne_get_shared_neuron_spktrain_subsets(exp_site_nedata, varargin)

% Function calculates shared neuron spktrain subsets.
%   exp_site_nedata: standard NE analysis structure array
%
%   threshalpha: sets threshold of cNE activity
%
%   neuopt: to be used when one does not want to compute all shared neurons
%   but just a subset. Default is an empty vector for when all shared
%   neurons are to be computed.
%
%   spkthresh: Threshold for cNE activity to binarize cNE activity
%
% Updated 9/3/17 by JS. Cleaned up code and output is now just one
% structure array.
% Updated 7/23/18 by JS. Cleaned up code and made it compatible with new
% 'ne_get_neuronal_ensemble_spktrain'. Added spkthresh option.

ip = inputParser;
addParameter(ip, 'threshalpha', 99.5, @(x) x >= 99 && x <= 99.9) %specifies how cNE activity is thresholded
addParameter(ip, 'neuopt', [], @(x) isempty(x) || isvector(x)) %specifies whether specific neurons are looked at
addParameter(ip, 'spkthresh', 100, @(x) x >= 0)
parse(ip, varargin{:})
threshalpha = ip.Results.threshalpha;
neuopt = ip.Results.neuopt;
spkthresh = ip.Results.spkthresh;


nedata = exp_site_nedata.nedata;
AM = cell2mat(nedata.NEmembers);

%get neurons which are members of more than 1 assembly
sorted = sort(AM);
nonuniqueidx = ~diff(sorted);
nonunique = sorted(nonuniqueidx);
nonunique = unique(nonunique);

if ~isempty(neuopt)
    %check if neurons selected are in more than one assembly
    Lia = ismember(neuopt, nonunique);
    wrongneu = neuopt(Lia == 0);
    if ~isempty(wrongneu)
        for ii = 1:length(wrongneu)
            warning(['Neuron %d does not belong to more than one '...
                'assembly, removing from output.'], wrongneu(ii))
        end
    end
    nonunique = neuopt(Lia);
        
end

if isempty(nonunique)
    warning('No shared neurons in this dataset');
    spkstruct = [];
    return
end

assemidx = cell(length(nonunique),1);
snlt = cell(size(assemidx)); 


NE_train = ne_get_neuronal_ensemble_spktrain(exp_site_nedata, 'threshalpha', threshalpha, 'method', 'repeat', 'memneuopt', 'none');

%convert spktrain to logical
if isfield(nedata, 'sta_spktrain')
    spktrain = nedata.sta_spktrain;
else
    spktrain = nedata.spktrain;
end
spklogic = spktrain>=1;

%get spktrains of neuron and assemblies that it is part of
for i = 1:length(nonunique)
    temp = cellfun(@(x) find(x == nonunique(i)),nedata.NEmembers,'UniformOutput',0);
    assemidx{i} = find(cellfun(@(x) (~isempty(x)),temp));
    snlt{i} = [spklogic(nonunique(i),:); NE_train(assemidx{i},:)];
end


spkstruct(length(nonunique)).neuron = [];

for j = 1:length(snlt)
   
    spkstruct(j).neuron = nonunique(j);
    spkstruct(j).NEs = assemidx{j};
    
    total = sum(snlt{j});
    
    % find when neuron spikes exclusively with each cNE
    withexc_idx = zeros(size(snlt{j},1)-1,length(spktrain));    
    for k = 2:size(snlt{j},1)
        withexc_idx(k-1,:) = (snlt{j}(k,:) == 1 & snlt{j}(1,:) == 1 & total == 2);
    end
            
    if spkthresh > 0
        
        tokeep = sum(withexc_idx, 2) >= spkthresh;
        if sum(tokeep) < 2 % remove neurons that are not part of more than one cNE after filtering            
            continue
        end
        
        % update struct array
        spkstruct(j).NEs_original = assemidx{j};
        spkstruct(j).NEs = assemidx{j}(tokeep);
        snlt{j} = snlt{j}(logical([1;tokeep]),:);
        
        % recalculate withexc_idx if any subsets were removed due to low
        % spike count
        if any(~tokeep)
            total = sum(snlt{j});             
            withexc_idx = zeros(size(snlt{j},1)-1,length(spktrain));    
            for k = 2:size(snlt{j},1)
                withexc_idx(k-1,:) = (snlt{j}(k,:) == 1 & snlt{j}(1,:) == 1 & total == 2);
            end
        end
        
    end
    
    % snlt: shared neurons logical train. Top row: neuronal logical spike train.
    % Subsequent rows: cNE-trains (repeated upsampling)
    spkstruct(j).snlt = snlt{j};
    
    %find when neuron spikes without cNE
    without_idx = total == 1 & snlt{j}(1,:) == 1;
    %find when neuron spikes with 1 or more cNEs
    with_idx = snlt{j}(1,:) == 1 & total >= 2;
    %find when neuron spikes with all cNEs
    all_idx = (total == size(snlt{j},1)); 
    
    spkstruct(j).spktrain = spktrain(nonunique(j),:);
    spkstruct(j).spktrain_count = sum(spkstruct(j).spktrain);
    spkstruct(j).withoutNE = spktrain(nonunique(j),:) .* without_idx;
    spkstruct(j).withoutNE_count = sum(spkstruct(j).withoutNE,2);
    spkstruct(j).withNE = spktrain(nonunique(j),:) .* with_idx;
    spkstruct(j).withNE_count = sum(spkstruct(j).withNE,2);
    spkstruct(j).withallNEs = spktrain(nonunique(j),:).* all_idx;
    spkstruct(j).withallNEs_count = sum(spkstruct(j).withallNEs,2);    
    spkstruct(j).NEexclusive = repmat(spktrain(nonunique(j),:),size(withexc_idx,1),1).* withexc_idx;
    spkstruct(j).NEexclusive_count = sum(spkstruct(j).NEexclusive,2);

end

% remove empty spkstructs
if spkthresh > 0
    spkstruct(cellfun('isempty', {spkstruct.snlt})) = [];
end
    
