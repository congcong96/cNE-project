function exp_site_nedata = ne_sort_exp_site_nedata_based_on_old(newdata, olddata)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

newpat = newdata.nedata.Patterns;
oldpat = olddata.nedata.Patterns;

corrmat = abs(corr(newpat, oldpat));
[~, idx] = max(corrmat, [], 1); 

exp_site_nedata = newdata;

exp_site_nedata.nedata.Patterns = exp_site_nedata.nedata.Patterns(:, idx);
exp_site_nedata.nedata.Activities = exp_site_nedata.nedata.Activities(idx, :);
exp_site_nedata.nedata.NEmembers = exp_site_nedata.nedata.NEmembers(idx);
exp_site_nedata.nedata.ca_stamat = exp_site_nedata.nedata.ca_stamat(idx,:);
exp_site_nedata.nedata.pwc = exp_site_nedata.nedata.pwc(idx);

end

