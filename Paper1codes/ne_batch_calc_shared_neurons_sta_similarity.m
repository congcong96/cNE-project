function ne_batch_calc_shared_neurons_sta_similarity(files)

for i = 1:length(files)
    
    load(files{i}, 'exp_site_nedata')
    clc;
    fprintf('\nProcessing %s...\n',files{i})
    
    shared_STAcorr = ne_calc_shared_neurons_sta_similarity(exp_site_nedata, 'sigopt', 0);
    save(files{i}, 'shared_STAcorr', '-append')
    clear('shared_STAcorr', 'exp_site_nedata')
       
end
