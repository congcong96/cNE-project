function p = ne_plot_ne_vs_neuron_info(files, infotype)

NEinfo = [];
neuroninfo = [];

for i = 1:length(files)
    s = load(files{i}, infotype);
    
    NEinfo = [NEinfo s.(sprintf('%s',infotype)).NE_info_extrap];
    try
        neuroninfo = [neuroninfo s.(sprintf('%s',infotype)).neuron_info_extrap];
    catch
        neuroninfo = [neuroninfo s.(sprintf('%s',infotype)).pNE_info_extrap];
    end

    
end

figure; hold on;
scatter(neuroninfo, NEinfo, 20, '.');
h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
xlabel('neuron MI (bits/spike)')
ylabel('NE MI (bits/event)')

[~, p(1)] = ttest(neuroninfo, NEinfo);
p(2) = signrank(neuroninfo, NEinfo);

return


    