function p = ne_plot_NE_vs_pNE_info(files)

NEinfo = [];
pNEinfo = [];

for i = 1:length(files)
    load(files{i}, 'pseudogroupinfo')
    
    NEinfo = [NEinfo pseudogroupinfo.NE_info_extrap];
    pNEinfo = [pNEinfo pseudogroupinfo.pNE_info_extrap];
    
end

figure; hold on;
scatter(pNEinfo, NEinfo);
h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
xlabel('neuron MI (bits/spike)')
ylabel('NE MI (bits/event)')

p(1) = ranksum(pNEinfo,NEinfo);
[~, p(2)] = ttest2(pNEinfo,NEinfo);

return


    