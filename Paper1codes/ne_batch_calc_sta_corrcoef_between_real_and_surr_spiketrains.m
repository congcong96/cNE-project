function [corrvals, staPred] = ne_batch_calc_sta_corrcoef_between_real_and_surr_spiketrains(files)

for i = 1:length(files)
    load(files{i})
    [corrvals{i}, staPred{i}] = ne_calc_sta_corrcoef_between_real_and_surr_spiketrains(exp_site_nedata, 0);
end