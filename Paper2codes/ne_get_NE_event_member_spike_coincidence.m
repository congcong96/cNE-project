function tab = ne_get_NE_event_member_spike_coincidence(exp_site_nedata)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;
NEmembers = nedata.NEmembers;

NEtrain = downsample_spiketrain(nedata.sta_NEtrain, 2);
spktrain = nedata.spktrain;
tab = cell(length(NEmembers), 1);

for i = 1:length(NEmembers)
    
    tempspktrain = spktrain(NEmembers{i},:);
    tempidx = tempspktrain(:, logical(NEtrain(i,:)));
    
    idxcell = cell(size(tempidx, 2), 1);
    
    for j = 1:size(tempidx, 2)
        
        idxcell{j} = find(tempidx(:,j));
    end
    
    idxcell = cellfun(@(x) NEmembers{i}(x), idxcell, 'UniformOutput', 0);
    [uniq, ~, idx] = unique(cellfun(@mat2str, idxcell, 'UniformOutput', 0));
    
    uniq = strrep(uniq, 'zeros(0,1)', '0');
    count = rude(sort(idx));
    
    tab{i} = table(count(:),'RowNames', uniq(:), 'VariableNames', {'Count'});
    if ~any(strcmp(tab{i}.Properties.RowNames, '0'))
        tab{i}{'0','Count'} = 0;
    end
end

