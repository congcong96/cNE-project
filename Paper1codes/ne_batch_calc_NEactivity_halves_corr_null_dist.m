function [actcorrcell, shuffledactcorrcell] = ne_batch_calc_NEactivity_halves_corr_null_dist(files)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

actcorrcell = cell(length(files),1);
shuffledactcorrcell = cell(length(files),1);

for i = 1:length(files)
    clc;
    fprintf('Processing %s...\n', files{i})
    load(files{i}, 'exp_site_nedata')
    
    [actcorrcell{i}, shuffledactcorrcell{i}] = ne_calc_NEactivity_halves_corr_null_dist(exp_site_nedata);


end

