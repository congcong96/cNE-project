function resultstruct = ne_calc_NE_statistics_by_group(sig_sta, statopt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

destructive_idx = cellfun(@(x) strcmp(x, 'destructive'), {sig_sta.category});
stimindependent_idx = cellfun(@(x) strcmp(x, 'stimindependent'), {sig_sta.category});
constructive_idx = cellfun(@(x) strcmp(x, 'constructive'), {sig_sta.category});
facilitative_idx = cellfun(@(x) strcmp(x, 'facilitative'), {sig_sta.category});

groups.facilitative = sig_sta(facilitative_idx);
groups.constructive = sig_sta(constructive_idx);
groups.independent = sig_sta(stimindependent_idx);
groups.destructive = sig_sta(destructive_idx);

fn = fieldnames(groups);

for i = 1:length(fn)
    
    tempstruct = groups.(fn{i});
    uniquefiles = unique({tempstruct.filename});
    tempcell = cell(length(uniquefiles), 1);
    
    for j = 1:length(uniquefiles)
        
        load(uniquefiles{j}, 'exp_site_nedata')
        idx = strcmp(uniquefiles{j}, {tempstruct.filename});
        switch statopt
            case 'depth'
                neurons = cell2mat({tempstruct(idx).neurons}');
                tempcell{j} = cellfun(@(x) x(2), exp_site_nedata.nedata.position(neurons)');
            case 'pairwise_distance'
                neurons = {tempstruct(idx).neurons}';
                temptempcell = cell(length(neurons), 1);
                for k = 1:length(neurons)
                    comb = nchoosek(neurons{k}, 2);
                    temptempcell{k} = ne_calc_interneuronal_distances(comb,...
                        cell2mat(exp_site_nedata.nedata.position'), 'distance');
                end
                tempcell{j} = cell2mat(temptempcell);
            case 'span'
                neurons = {tempstruct(idx).neurons}';
                pos = cellfun(@(x) x(2), exp_site_nedata.nedata.position);
                depth = cellfun(@(x) pos(x), neurons, 'UniformOutput', 0);
                tempcell{j} = cellfun(@(x) max(x) - min(x), depth);
            case 'median'
                neurons = {tempstruct(idx).neurons}';
                pos = cellfun(@(x) x(2), exp_site_nedata.nedata.position);
                depth = cellfun(@(x) pos(x), neurons, 'UniformOutput', 0);
                tempcell{j} = cellfun(@median, depth);
            case 'percentage_coincidence'
                NEs = [tempstruct(idx).NE]';
                NEmemtrain = ne_get_neuronal_ensemble_spktrain(exp_site_nedata,...
                    'threshalpha', 99.5, 'method', 'repeat', 'memneuopt', 'raw');
                NEmemtrain = NEmemtrain(NEs);
                tempcell{j} = cell2mat(cellfun(@(x) sum(x(2:end, logical(x(1,:))), 2)...
                    ./ sum(x(2:end, :), 2), NEmemtrain, 'UniformOutput',0));
            case 'percentage_coincidence_cNE_activity'
                NEs = [tempstruct(idx).NE]';
                NEtrain = downsample_spiketrain(exp_site_nedata.nedata.sta_NEtrain, 2);
                NEtrain = NEtrain(NEs,:);
                tempcell{j} = zeros(size(NEtrain, 1), 1);
                for k = 1:size(NEtrain, 1)
                    combinedtrain = [NEtrain(k,:); exp_site_nedata.nedata.spktrain(exp_site_nedata.nedata.NEmembers{NEs(k)}, :)];
                    tempcell{j}(k) = sum(combinedtrain(1,:) & sum(combinedtrain(2:end,:), 1) > 0) ./ ...
                        sum(combinedtrain(1,:));
                end
        end
    end
    resultstruct.(fn{i}) = cell2mat(tempcell);
    
end    

end

