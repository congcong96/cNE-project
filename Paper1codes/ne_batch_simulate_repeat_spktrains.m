function [realspktrains, surrspktrains, varargout] = ne_batch_simulate_repeat_spktrains(files, niter, goodneurons)

% Simulates spiketrains of n neurons where n is number of repeats recorded

if nargin == 2
    goodneurons = [];
end

spktrains = cell(length(files),1);

% Get all spiketrains
for i = 1:length(files)
    load(files{i})
    spktrains{i} = exp_site_nedata.nedata.spktrain;    
end

realspktrains = cell(niter,1);
surrspktrains = cell(niter,1);

for j = 1:niter
    if isempty(goodneurons)
        comb = randsample(size(spktrains{1},1), length(files));
    else
        comb = randsample(goodneurons, length(files));
    end
    varargout{2}(j,:) = comb;
    realspktrains{j} = cell(length(files),1);
    surrspktrains{j} = cell(length(files),1);
    
    for k = 1:length(files)
        [realcomb, surrcombidx] = sort(comb); %sort neuronal identity in numerical order
        realspktrains{j}{k} = spktrains{k}(realcomb,:); 
        surrspktrains{j}{k} = cell2mat(arrayfun(@(x,y) spktrains{x}(y,:),1:length(files),...
            comb','UniformOutput', 0)');
        surrspktrains{j}{k} = surrspktrains{j}{k}(surrcombidx,:); %sort neuronal identity in numerical order
        
%         %test
%         sumtest = sum(surrspktrains{j}{k},2);
%         for ii = 1:length(comb)
%             assert(sum(spktrains{ii}(comb(ii),:)) == sumtest(ii == surrcombidx));
%         end
        
        comb = circshift(comb,1);
    end
    varargout{1}(j,:) = realcomb;
end
        

end

