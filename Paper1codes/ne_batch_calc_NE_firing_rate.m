function rate = ne_batch_calc_NE_firing_rate(NEfiles, alpha)


if ~exist('alpha','var')
    alpha = [90 92.5 95 97.5 99 99.5 99.9 99.99];
end

rate = cell(1, length(NEfiles));

for i = 1:length(NEfiles)
    clc; fprintf('Processing %s\n\n', NEfiles{i});
    load(NEfiles{i}, 'exp_site_nedata')
    [rate{i}] = ne_calc_NE_firing_rate(exp_site_nedata, alpha);
   
end

