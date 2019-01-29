function [NEcoin, PSTH] = ne_calc_MU_coincidence_with_NE_events(exp_site_nedata, MUthresh, trigger, stim_mat)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('stim_mat','var')
    stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata, 1);
    spktrain = ne_create_spktrain_from_stim_mat(MUthresh, stimstr.stimulus, trigger);
else
    spktrain = ne_create_spktrain_from_stim_mat(MUthresh, stim_mat, trigger);
end

spktrain = downsample_spiketrain(spktrain, exp_site_nedata.df);
PSTH = zscore(sum(spktrain, 1));

nedata = exp_site_nedata.nedata;
NEact = nedata.Activities;

alpha = 99.5;
NEthresh = nedata.NEthresh((nedata.NEthresh_alpha == alpha), :);

NEcoin = cell(size(NEact,1), 1);
for i = 1:size(NEact, 1)
    NEtrain = NEact(i,:) >= NEthresh(i);
    NEcoin{i} = PSTH(NEtrain);% ./ sum(PSTH); %length(PSTH));
end

