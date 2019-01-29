function ne_batch_group_classified_STAs_based_on_NE(nefiles)
%UNTITLED24 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(nefiles)
    
    fprintf('\nGrouping STA siginficance for %s...\n', nefiles{i})
    
    load(nefiles{i}, 'exp_site_nedata')
    sig_sta = ne_group_classified_STAs_based_on_NE(exp_site_nedata);
    
    exp_site_nedata.nedata.sig_sta = sig_sta;

    save(nefiles{i}, 'exp_site_nedata', '-append');
    clear('exp_site_nedata')
    
end

end

