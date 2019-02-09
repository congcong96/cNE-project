function [corrcell, shuffledcorrcell] = ne_batch_calc_ICweight_halves_null_dist(files)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

corrcell = cell(length(files),1);
shuffledcorrcell = cell(length(files),1);

for i = 1:length(files)
    
    clc;
    fprintf('Processing %s...\n', files{i})
    
    load(files{i}, 'exp_site_nedata')
    [corrcell{i}, shuffledcorrcell{i}] = ne_calc_IC_halves_null_dist(exp_site_nedata);


end

