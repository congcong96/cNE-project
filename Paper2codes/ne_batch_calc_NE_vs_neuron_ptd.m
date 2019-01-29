function ne_batch_calc_NE_vs_neuron_ptd(nefiles)

for i = 1:length(nefiles)
    
    clc;
    fprintf('Processing files %d of %d...\n', i, length(nefiles))
    
    load(nefiles{i}, 'exp_site_nedata')
    NEneuronptdall = ne_calc_NE_vs_neuron_ptd(exp_site_nedata);
    
    save(nefiles{i}, 'NEneuronptdall', '-append');
end