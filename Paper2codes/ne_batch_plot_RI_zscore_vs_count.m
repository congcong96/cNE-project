function ne_batch_plot_RI_zscore_vs_count(nefiles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

NE_RI_zscore = cell(length(nefiles), 1);
neuron_RI_zscore = cell(length(nefiles), 1);
NE_count = cell(length(nefiles), 1);
neuron_count = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata');   
    NE_RI_zscore{i} = exp_site_nedata.nedata.sig_classifying_stats.NE_reliability_idx_zscore;
    neuron_RI_zscore{i} = exp_site_nedata.nedata.sig_classifying_stats.neuron_reliability_idx_zscore;
    NE_count{i} = sum(exp_site_nedata.nedata.sta_NEtrain, 2);
    neuron_count{i} = sum(exp_site_nedata.nedata.sta_spktrain, 2);

end

figure;
scatter(cell2mat(NE_count), cell2mat(NE_RI_zscore), 10);

figure;
scatter(cell2mat(neuron_count), cell2mat(neuron_RI_zscore), 10);