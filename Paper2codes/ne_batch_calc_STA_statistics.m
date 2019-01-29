function ne_batch_calc_STA_statistics(nefiles)

for i = 1:length(nefiles)
    
    fprintf('\nProcessing STA statistics for %s...\n', nefiles{i})
    load(nefiles{i}, 'exp_site_nedata')
    STA_stats = ne_calc_STA_statistics_from_exp_site_nedata(exp_site_nedata);
    exp_site_nedata.nedata.STA_stats = STA_stats;
    
    save(nefiles{i}, 'exp_site_nedata', '-append');
end
    
