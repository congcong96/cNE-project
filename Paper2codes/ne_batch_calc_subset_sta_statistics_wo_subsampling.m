function ne_batch_calc_subset_sta_statistics_wo_subsampling(nefiles)

for i = 1:length(nefiles)
    clc;
    fprintf('Processing file %d of %d...\n', i, length(nefiles))
    
    load(nefiles{i}, 'exp_site_nedata', 'spksubset')
    nosample_subsetstats = ne_calc_subset_sta_statistics_wo_subsampling(exp_site_nedata, spksubset);
    save(nefiles{i}, 'nosample_subsetstats','-append');

end

