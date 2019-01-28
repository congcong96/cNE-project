function ne_batch_calc_NE_act_thresholds(files)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

alpha = 99:0.1:99.9;

for i = 1:length(files)
    load(files{i})
    thresh = ne_calc_NE_act_thresholds(exp_site_nedata,'circular', 50, alpha);
    exp_site_nedata.nedata.NEthresh = thresh;
    exp_site_nedata.nedata.NEthresh_alpha = alpha;
    save(files{i}, 'exp_site_nedata', '-append');

end

