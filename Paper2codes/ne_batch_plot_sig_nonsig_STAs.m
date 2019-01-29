function ne_batch_plot_sig_nonsig_STAs(nefiles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(nefiles)
    load(nefiles{i})
    ne_plot_sig_nonsig_STAs(exp_site_nedata);
    
    pause;
    close all


end

