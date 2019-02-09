function [FSI] = ne_batch_calc_NE_neuron_FSI(files)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

FSI_ne = cell(length(files),1);
FSI_neuron = cell(length(files),1);
FSI = cell(length(files),1);
for i = 1:length(files)
    load(files{i}, 'exp_site_nedata')
    fprintf('\nCalculating FSI for %s...', files{i})
%     [FSI_ne{i}, FSI_neuron{i}] = ne_calc_NE_neuron_FSI(exp_site_nedata);
    FSI{i} = ne_calc_NE_neuron_paired_FSI(exp_site_nedata);

end

fprintf('\n')
