function p = ne_batch_plot_sig_vs_nonsig_neuron_depth(nefiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sigdepth = cell(length(nefiles), 1);
nonsigdepth = cell(length(nefiles), 1);
color = eight_color_blind_palette('orange', 'blue');


for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    
    sigdepth{i} = cellfun(@(x) x(2), exp_site_nedata.nedata.position(...
        exp_site_nedata.nedata.sig_neuron_sta));
    
    nonsigdepth{i} = cellfun(@(x) x(2), exp_site_nedata.nedata.position(...
        ~exp_site_nedata.nedata.sig_neuron_sta));


end


depthstruct.sig = cell2mat(sigdepth')';
depthstruct.nonsig = cell2mat(nonsigdepth')';

figure('Position', [680 478 372 500]);
axis ij
p = plot_plotspread_and_boxplot(depthstruct, 'ranksum', color, -1);
ylabel('Neuronal depth (um)')
print_mfilename(mfilename);