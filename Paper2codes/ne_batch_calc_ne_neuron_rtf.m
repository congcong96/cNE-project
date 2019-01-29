function ne_batch_calc_ne_neuron_rtf(nefiles, sigopt)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(nefiles)
    
    fprintf('\nProcessing RTFs for %s...\n', nefiles{i})
    
    load(nefiles{i}, 'exp_site_nedata')
    [tmf, xmf, neuronrtf, NErtf] = ne_calc_ne_neuron_rtf(exp_site_nedata, sigopt);
     
     exp_site_nedata.nedata.tmf = tmf;
     exp_site_nedata.nedata.xmf = xmf;
     exp_site_nedata.nedata.neuron_rtfmat = neuronrtf;
     exp_site_nedata.nedata.NE_rtfmat = NErtf;
     
     save(nefiles{i}, 'exp_site_nedata', '-append');
     clear('exp_site_nedata')
     
end

