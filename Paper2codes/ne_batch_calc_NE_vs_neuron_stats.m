function ne_batch_calc_NE_vs_neuron_stats(nefiles)

for i = 1:length(nefiles)
    
    clc;
    fprintf('Processing files %d of %d...\n', i, length(nefiles))
    
    load(nefiles{i}, 'exp_site_nedata')
    NEneuronstats = ne_calc_NE_vs_neuron_stats(exp_site_nedata);
    
    save(nefiles{i}, 'NEneuronstats', '-append');
end