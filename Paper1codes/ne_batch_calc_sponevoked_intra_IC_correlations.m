function [sponcell, evokedcell] = ne_batch_calc_sponevoked_intra_IC_correlations(files)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sponcell = cell(length(files),1);
evokedcell = cell(length(files),1);

for i = 1:length(files)
    
    load(files{i}, 'exp_site_nedata')
    [sponcell{i}, evokedcell{i}] = ne_calc_sponevoked_intra_IC_correlations(exp_site_nedata);
    

end

