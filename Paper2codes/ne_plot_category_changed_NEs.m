function ne_plot_category_changed_NEs(sig_sta1, sig_sta2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

category1 = {sig_sta1.category};
category2 = {sig_sta2.category};

idx = find(~cellfun(@(x,y) strcmp(x,y), category1, category2));


changedcat1 = category1(idx);
changedcat2 = category2(idx);

changedNE = sig_sta1(idx);


for i = 1:length(changedNE)
    
     ne_plot_ensemble_members_sta_from_sigsta(changedNE(i), 1);
     suptitle(sprintf('index:%d, %s -> %s', idx(i), changedcat1{i}, changedcat2{i}))

end

