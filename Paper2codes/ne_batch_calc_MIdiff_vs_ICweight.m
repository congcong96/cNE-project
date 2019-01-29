function [bcell, pvalcell] = ne_batch_calc_MIdiff_vs_ICweight(nefiles)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

bcell = cell(length(nefiles), 1);
pvalcell = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata', 'NEneuroninfo')
    
    if isempty(NEneuroninfo)
        continue
    end
    
    [bcell{i}, pvalcell{i}] = ne_plot_MIdiff_vs_ICweight(exp_site_nedata, NEneuroninfo, 0);


end

