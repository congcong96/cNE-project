function statstruct = ne_calc_subset_sta_statistics(subsetstruct, subsampleopt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

stims = regexp({subsetstruct.filename},'(?<=db-)rn\d{1,2}-\d{2}min(?=-[aH])','match','once');
uniq_stims = unique(stims);

fn = fieldnames(subsetstruct);
spkfn = fn(~cellfun('isempty', regexp(fn, 'spktrain$')));

if subsampleopt
    noncountfn = fn(cellfun('isempty', regexp(fn, 'count$')));
    countarray = squeeze(cell2mat(struct2cell(rmfield(subsetstruct, noncountfn))));
    minspikes = floor(0.99*min(countarray,[],1));
end

statstruct = rmfield(subsetstruct, spkfn);

for i = 1:length(uniq_stims)
    
    idx = find(contains(stims, uniq_stims{i}));    
    load(subsetstruct(idx(1)).filename, 'exp_site_nedata');
    stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);
    
    if subsampleopt
        tempmin = minspikes(idx);
    end        
    
    for ii = 1:length(spkfn)
        
        spktemp = cell2mat({subsetstruct(idx).(spkfn{ii})}');
        prefix = regexp(spkfn{ii}, '^\w+(?=_)','match','once');
    
        if subsampleopt
            niter = 20;
            spktrain_subsample = cell(niter, 1);
            
            for j = 1:niter            
                spktrain_subsample{j} = zeros(size(spktemp));            
                for k = 1:size(spktemp, 1)                
                    spktrain_subsample{j}(k,:) = sub_sample_spktrain(spktemp(k,:), sum(spktemp(k,:)) - tempmin(k));
                end
            end
        
            spktrain_subsample = cell2mat(spktrain_subsample);
            stastats = ne_calc_sta_stats_zscore_from_spktrain(spktrain_subsample, stimstr.stimulus, 5);
            fnstats = fieldnames(stastats);
            newfn = cellfun(@(x) [prefix '_' x], fnstats, 'UniformOutput', 0);       

            for j = 1:length(fnstats)
                temp = stastats.(fnstats{j});
                for k = 1:size(spktemp, 1)
                    temptemp = temp(k:size(spktemp,1):end, :);
                    statstruct(idx(k)).(newfn{j}) = temptemp;
                end
            end

        else        
            stastats = ne_calc_sta_stats_zscore_from_spktrain(spktemp, stimstr.stimulus);
            fnstats = fieldnames(stastats);
            newfn = cellfun(@(x) [prefix '_' x], fnstats, 'UniformOutput', 0);        


            for j = 1:length(fnstats)
                temp = mat2cell(stastats.(fnstats{j}), ones(size(stastats.(fnstats{j}),1), 1), size(stastats.(fnstats{j}), 2));
                [statstruct(idx).(newfn{j})] = temp{:};
            end

        end      
    
    end
end


