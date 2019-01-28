function ne_batch_sort_exp_site_nedata_based_on_old(oldfiles, newfiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

assert(length(oldfiles) == length(newfiles))

for i = 1:length(oldfiles)
    load(oldfiles{i})
    olddata = exp_site_nedata;
    clear('exp_site_nedata')
    
    load(newfiles{i})
    newdata = exp_site_nedata;
    clear('exp_site_nedata')
    
    exp_site_nedata = ne_sort_exp_site_nedata_based_on_old(newdata, olddata);
    save(newfiles{i}, 'exp_site_nedata')

end

