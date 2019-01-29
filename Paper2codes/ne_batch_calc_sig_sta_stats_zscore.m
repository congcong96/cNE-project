function ne_batch_calc_sig_sta_stats_zscore(nefiles)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(nefiles)
    
    fprintf('\nProcessing %s (files %d of %d)...',nefiles{i},i,length(nefiles))
    load(nefiles{i})
    nedata = exp_site_nedata.nedata;
    
    STA_stats = nedata.STA_stats;
    shuffled_stats = nedata.shuffled_stats;
    
    fn = fieldnames(nedata.STA_stats);
    
    for j = 1:length(fn)
                
        miu = mean(shuffled_stats.(fn{j}), 2, 'omitnan');
        sigma = std(shuffled_stats.(fn{j}),0, 2, 'omitnan');
        nedata.sig_classifying_stats.([fn{j} '_zscore']) = (mean(STA_stats.(fn{j}), 2, 'omitnan') - miu) ./ sigma;
        
    end
    
    exp_site_nedata.nedata = nedata;
    save(nefiles{i}, 'exp_site_nedata', '-append');
    clear('exp_site_nedata')
    

end

fprintf('\n')