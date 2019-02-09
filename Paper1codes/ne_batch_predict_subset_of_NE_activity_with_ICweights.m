function corractmat = ne_batch_predict_subset_of_NE_activity_with_ICweights(files)

corractmat = cell(length(files),1);

for i = 1:length(files)
    
    load(files{i}, 'exp_site_nedata')
    
    corractmat{i} = ne_predict_subset_of_NE_activity_with_ICweights(exp_site_nedata);

end


