function statstruct = ne_calc_subset_sta_statistics_wo_subsampling(exp_site_nedata, spksubset)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% stims = regexp({subsetstruct.filename},'(?<=db-)rn\d{1,2}-\d{2}min(?=-[aH])','match','once');
% uniq_stims = unique(stims);
% 
% statstruct = subsetstruct;
statstruct = spksubset;

for i = 1:length(uniq_stims)
    
    idx = find(contains(stims, uniq_stims{i}));    
    load(spksubset(idx(1)).filename, 'exp_site_nedata');
    stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);
    
    spktrain = cell2mat({spksubset(idx).spktrain}');
    subset_spktrain = cell2mat({spksubset(idx).subset_spktrain}');
    
    if subsampleopt
        niter = 20;
        subset_spkcount = sum(subset_spktrain, 2);
        spktrain_spkcount = sum(spktrain, 2);
        spktrain_subsample = cell(niter, 1);
        
        for j = 1:niter            
            spktrain_subsample{j} = zeros(size(spktrain));            
            for k = 1:size(spktrain, 1)                
                spktrain_subsample{j}(k,:) = sub_sample_spktrain(spktrain(k,:), spktrain_spkcount(k) - subset_spkcount(k));
            end
        end
        
        spktrain_subsample = cell2mat(spktrain_subsample);
        stastats = ne_calc_sta_stats_zscore_from_spktrain(spktrain_subsample, stimstr.stimulus);
        fn = fieldnames(stastats);
        newfn = cellfun(@(x) ['all_' x], fn, 'UniformOutput', 0);  
        
        for j = 1:length(fn)
            temp = stastats.(fn{j});
            for k = 1:size(spktrain, 1)
                temptemp = temp(k:size(spktrain,1):end, :);
                statstruct(idx(k)).(newfn{j}) = temptemp;
            end
        end
        
    else        
        stastats = ne_calc_sta_stats_zscore_from_spktrain(spktrain, stimstr.stimulus);
        fn = fieldnames(stastats);
        newfn = cellfun(@(x) ['all_' x], fn, 'UniformOutput', 0);        
    

        for j = 1:length(fn)
            temp = mat2cell(stastats.(fn{j}), ones(size(stastats.(fn{j}),1), 1), size(stastats.(fn{j}), 2));
            [statstruct(idx).(newfn{j})] = temp{:};
        end
        
    end
    
    stastats = ne_calc_sta_stats_zscore_from_spktrain(subset_spktrain, stimstr.stimulus);
    fn = fieldnames(stastats);
    newfn = cellfun(@(x) ['withNE_' x], fn, 'UniformOutput', 0);
    
    for j = 1:length(fn)
        temp = mat2cell(stastats.(fn{j}), ones(size(stastats.(fn{j}),1), 1), size(stastats.(fn{j}), 2));
        [statstruct(idx).(newfn{j})] = temp{:};
    end
    
end

