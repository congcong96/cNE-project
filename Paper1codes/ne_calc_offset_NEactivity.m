function usNE = ne_calc_offset_NEactivity(spkfile, df, reqdf, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addOptional(p, 'NEthreshopt', 0, @(x) x == 0 | x == 1);
addOptional(p, 'stim_mat', [])
parse(p, varargin{:})
NEthreshopt = p.Results.NEthreshopt;
stim_mat = p.Results.stim_mat;

basebin = 0.5; % 1 dft = 0.5 ms time bins
n = df/reqdf; % number of times to calculate offset
offset = 0;

% load spikes and trigger and get spike times
load(spkfile, 'spk', 'trigger')
orispktimes = {spk.spiketimes};

% load stimulus matrix if required
if isempty(stim_mat)
    stim = regexp(spkfile,'(?<=(db-))rn\d{1,2}','match','once');

    drive = gcdr;
    subfolder = 'Ripple_Noise\downsampled';
    stimfolder = fullfile(drive,subfolder);

    stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stim);
    stimfile = fullfile(stimfolder,stimfilename);
    stimfile = dir(stimfile);
    load(stimfile.name);
end

% initialize struct array
usNE(n).offset = [];

% get NE for all offset datasets
spktimes = orispktimes;

for i = 1:n    
    
    % save stuff
    usNE(i).offset = offset;
    
    
    [spktrain, ~] = ca_create_spktrain_from_stim_mat(spktimes, stim_mat, trigger);
    spkmat = downsample_spiketrain(spktrain, df);
    [usNE(i).unsorted_pat, usNE(i).unsorted_act] = ca_detect_cell_assemblies_data(spkmat);
    
%     if i == 1
%         offset = basebin * reqdf;
%     elseif mod(i,2) == 0
%         offset = -offset;
%     else %mod(i,2) == 1
%         offset = abs(offset) + basebin * reqdf;
%     end
    offset = offset + basebin * reqdf;
    spktimes = cellfun(@(x) x + offset, orispktimes, 'UniformOutput', 0);
    
    close all
    
    usNE(i).spktrain = spkmat;

    
end

    
minNE = min(cellfun(@(x) size(x, 2), {usNE.unsorted_pat}));

% set first (and reference) NE to have the min number of NEs
if size(usNE(1).unsorted_pat,2) ~= minNE
    firstminNEidx = find(cellfun(@(x) size(x,2), {usNE.unsorted_pat}) == minNE, 1);
    tempcorrmat = abs(corr(usNE(1).unsorted_pat, usNE(firstminNEidx).unsorted_pat));
    [row, ~] = find(tempcorrmat > 0.75);
    usNE(1).Patterns = usNE(1).unsorted_pat(:,row);
    usNE(1).Activities = usNE(1).unsorted_act(row, :);
else
    usNE(1).Patterns = usNE(1).unsorted_pat;
    usNE(1).Activities = usNE(1).unsorted_act;
end


for j = 2:n

    corrmat = abs(corr(usNE(j).unsorted_pat, usNE(1).Patterns));
    [row, col] = find(corrmat > 0.75);
    
    usNE(j).Patterns = usNE(j).unsorted_pat(:,row);
    usNE(j).Activities = usNE(j).unsorted_act(row,:);
    
    usNE(j).corrmat = corr(usNE(j).Patterns, usNE(1).Patterns);
    
    
    if length(row) ~= minNE     
        usNE(j).unmatched_NE = setdiff(1:size(usNE(1).Patterns, 2), col);
    else
        usNE(j).unmatched_NE = [];
    end
        

end
  
% correct for unmatched NEs
unmatchedcell = {usNE.unmatched_NE};
unmatchedidx = ~cellfun('isempty', unmatchedcell);

if any(unmatchedidx)
    unmatched = unique(cell2mat(unmatchedcell(unmatchedidx)));
    
    for i = 1:length(usNE)
        if ~unmatchedidx(i)            
            usNE(i).Patterns(:,unmatched) = [];
            usNE(i).Activities(unmatched,:) = [];            
        end
        usNE(i).corrmat = corr(usNE(i).Patterns, usNE(1).Patterns);
    end
end

% set baseline activity to 0
for j = 1:length(usNE)
    for k = 1:size(usNE(j).Activities, 1)
        usNE(j).Activities(k,:) = usNE(j).Activities(k,:) - usNE(j).Activities(k,1);
    end
end

if NEthreshopt == 1
    for k = 1:length(usNE)
        usNE(k).NEthresh = ne_calc_NE_act_thresholds(usNE(k),'circular',50,99.9);
    end
end

   
end

