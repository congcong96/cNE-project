function [rate, alpha] = ne_calc_NE_firing_rate(exp_site_nedata, alpha)

if ~exist('alpha','var')
    alpha = [90 95 97.5 99 99.5 99.9 99.99];
end

df = exp_site_nedata.df;
nedata = exp_site_nedata.nedata;
act = nedata.Activities;
len = size(act, 2);

tottime = 0.5 * df * len / 1000; % in seconds

thresh = ne_calc_NE_act_thresholds(exp_site_nedata,'circular',50,alpha);

rate = zeros(size(thresh));

for i = 1:size(thresh,1)
    
    for j = 1:size(thresh,2)
        
        rate(i,j) = sum(act(j,:) > thresh(i,j)) ./ tottime;
        
    end
end
