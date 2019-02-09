function repNE = ne_calc_subsets_repeats_NE_stats(files, varargin)

% Get statistics of NE from subsets of spike trains to match the number of
% repeats, i.e. for 15 repeats of a stimulus, we draw 15 neurons from one
% stimulus (for each of the presentations) and compare against drawing 15
% neurons, one from each presentation. Neuronal identity must be preserved 
% within the whole dataset.

p = inputParser;
addRequired(p, 'files', @iscell);
addParameter(p, 'ICthresh',[], @isscalar); 
addParameter(p, 'goodneurons', [], @isvector);
parse(p, files, varargin{:});
files = p.Results.files;
ICthresh = p.Results.ICthresh;
goodneurons = p.Results.goodneurons;

[realspktrains, surrspktrains, neurons, neuronorder] = ...
    ne_batch_simulate_repeat_spktrains(files, 1, goodneurons);


% for i = 1:niter
%     clc;
%     fprintf('\nIteration %d of %d\n',i,niter)
L = length(files);

rICweights = cell(L,1);
sICweights = cell(L,1);
rCI = zeros(L,2);
sCI = zeros(L,2);
rnumNEs = zeros(L,1);
snumNEs = zeros(L,1);
rNEsize = zeros(L,1);
sNEsize = zeros(L,1);
rNEmembers = cell(L,1);
sNEmembers = cell(L,1);

for j = 1:L
    rICweights{j} = assembly_patterns(realspktrains{1}{j});
    sICweights{j} = assembly_patterns(surrspktrains{1}{j});
    if isempty(ICthresh)
        rCI(j,:) = ne_calc_ICA_threshold(realspktrains{1}{j}, 'circular', 50, 'stdev', 1.5);
        sCI(j,:) = ne_calc_ICA_threshold(surrspktrains{1}{j}, 'circular', 50, 'stdev', 1.5);
    else
        rCI(j,:) = [-ICthresh{1} ICthresh{1}];
        sCI(j,:) = [-ICthresh{1} ICthresh{1}];

    end

    rNEmembers{j} = cell(size(rICweights{j},2),1);
    sNEmembers{j} = cell(size(sICweights{j},2),1);

    rnumNEs(j) = size(rICweights{j},2);
    snumNEs(j) = size(sICweights{j},2);

    close all

    for ii = 1:size(rICweights{j},2)

        idx1 = find(rICweights{j}(:,ii) <= rCI(j,1));
        idx2 = find(rICweights{j}(:,ii) >= rCI(j,2));
        rNEmembers{j}{ii} = sort([idx1; idx2]);

    end

    for jj = 1:size(sICweights{j},2)

        idx3 = find(sICweights{j}(:,jj) <= sCI(j,1));
        idx4 = find(sICweights{j}(:,jj) >= sCI(j,2));
        sNEmembers{j}{jj} = sort([idx3; idx4]);

    end

    rNEsize(j) = mean(cellfun('length', rNEmembers{j}));
    sNEsize(j) = mean(cellfun('length', sNEmembers{j}));

end

repNE.neurons = neurons(1,:);
repNE.order = neuronorder(1,:);
repNE.realspktrains = realspktrains{1};
repNE.surrspktrains = surrspktrains{1};
repNE.realICwts = rICweights;
repNE.surrICwts = sICweights;
repNE.realCI = rCI;
repNE.surrCI = sCI;
repNE.realnumNEs = rnumNEs;
repNE.surrnumNEs = snumNEs;

repNE.realNEmembers = rNEmembers;
repNE.surrNEmembers = sNEmembers;
repNE.realNEsize = rNEsize;
repNE.surrNEsize = sNEsize;
    
%     save(sprintf('repNE%d',1), 'repNE')
%     
%     clear('repNE')
%     
% end
return
