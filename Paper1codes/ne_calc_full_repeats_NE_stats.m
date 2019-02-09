function ne_calc_full_repeats_NE_stats(files, niter)

spktrains = cell(length(files),1);

% Get all spiketrains
for i = 1:length(files)
    load(files{i})
    spktrains{i} = exp_site_nedata.nedata.spktrain;    
end

% Get number of neurons and repeats
numneurons = size(exp_site_nedata.nedata.spktrain,1);
numfiles = length(files);
clear('exp_site_nedata');

rem = mod(numneurons, numfiles);
reps = floor(numneurons/numfiles);

for j = 1:niter
    clc;
    fprintf('Running iteration %d...\n',j)
    
    extras = randsample(numfiles, rem);
    NEidx = [repmat(1:numfiles,1,reps)'; extras];
    NEidx = NEidx(randperm(numneurons));
    spktrain = zeros(size(spktrains{1}));
    
    for k = 1:numfiles
        neuidx = find(NEidx == k);
        spktrain(neuidx,:) = spktrains{k}(neuidx,:);
    end
    
    ICweights = assembly_patterns(spktrain);
    CI = ne_calc_ICA_threshold(spktrain, 'circular');
    
    NEmembers = cell(size(ICweights,2),1);

    for ii = 1:size(ICweights,2)

        idx1 = find(ICweights(:,ii) <= CI(1));
        idx2 = find(ICweights(:,ii) >= CI(2));
        NEmembers{ii} = sort([idx1; idx2]);

    end
    
    repNE.NEidx = NEidx;
    repNE.spktrain = spktrain;
    repNE.ICwts = ICweights;
    repNE.numNEs = size(ICweights,2);
    repNE.CI = CI;
    repNE.NEmembers = NEmembers;
    repNE.NEsize = mean(cellfun('length', NEmembers));
    
    save(sprintf('repNE%d',j), 'repNE')
    
    clear('repNE')
    close all
    
end

fprintf('\n')
        