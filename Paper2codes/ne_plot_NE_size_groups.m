function p = ne_plot_NE_size_groups(sig_sta)

destructive_idx = ([sig_sta.sig_NE] == 0 & [sig_sta.percentage_sig_neurons] > 0.5);
stimindependent_idx = ([sig_sta.sig_NE] == 0 & [sig_sta.percentage_sig_neurons] <= 0.5);
constructive_idx = ([sig_sta.sig_NE] == 1 & [sig_sta.percentage_sig_neurons] < 0.5);
facilitative_idx = ([sig_sta.sig_NE] == 1 & [sig_sta.percentage_sig_neurons] >= 0.5);

numneurons.facilitative = cellfun('length', {sig_sta(facilitative_idx).sig_neurons}');
numneurons.constructive = cellfun('length', {sig_sta(constructive_idx).sig_neurons}');
numneurons.independent = cellfun('length', {sig_sta(stimindependent_idx).sig_neurons}');
numneurons.destructive = cellfun('length', {sig_sta(destructive_idx).sig_neurons}');

color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');

figure;
p{1} = plot_plotspread_and_boxplot(numneurons, 'kruskalwallis', color, 1);
ylabel('Number of neurons in each cNE')

print_mfilename(mfilename);

percneurons.facilitative = cellfun('length', {sig_sta(facilitative_idx).sig_neurons}') ./...
    cellfun('length', {sig_sta(facilitative_idx).IC_weights}');
percneurons.constructive = cellfun('length', {sig_sta(constructive_idx).sig_neurons}') ./...
    cellfun('length', {sig_sta(constructive_idx).IC_weights}');
percneurons.independent = cellfun('length', {sig_sta(stimindependent_idx).sig_neurons}') ./...
    cellfun('length', {sig_sta(stimindependent_idx).IC_weights}');
percneurons.destructive = cellfun('length', {sig_sta(destructive_idx).sig_neurons}') ./...
    cellfun('length', {sig_sta(destructive_idx).IC_weights}');

figure;
p{2} = plot_plotspread_and_boxplot(percneurons, 'kruskalwallis', color, 1);
ylabel('Percentage of recorded neurons in each cNE')

print_mfilename(mfilename);

end

