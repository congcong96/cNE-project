function [coinratio, varargout] = ne_calc_coinratio_from_simulated_spktrains(files, model, comb, mtype)

allcoinratio = cell(length(comb), length(files));

for i = 1:length(files)
    fprintf('\nProcessing %s...\n',files{i})
    load(files{i})
    spiketrains = surrspktrain.(model);
    
    temp = cell(length(comb), 1);
    
    for j = 1:length(spiketrains)
        fprintf('%d of %d spiketrains\n', j, length(spiketrains))
        for k = 1:length(comb)
            
            temp{k}(:,j) = ne_calc_coincidence_within_spktrain(spiketrains{j}, comb{k});
            
        end           
        
    end
    
    allcoinratio(:,i) = temp;
        
end

coinratio = cell(size(allcoinratio,1),1);
for ii = 1:size(allcoinratio,1)
    switch mtype
        case 'median'
            coinratio{ii} = median(cell2mat(allcoinratio(ii,:)),2);
        case 'mean'
            coinratio{ii} = mean(cell2mat(allcoinratio(ii,:)),2);
    end
end

varargout{1} = allcoinratio;