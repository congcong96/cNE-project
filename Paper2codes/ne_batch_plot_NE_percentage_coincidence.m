function p = ne_batch_plot_NE_percentage_coincidence(sig_sta)

%% Box plot median depth all groups

medianstruct = ne_calc_NE_statistics_by_group(sig_sta, 'percentage_coincidence');

color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');

figure;
p = plot_plotspread_and_boxplot(medianstruct, 'kruskalwallis', color, 1);
print_mfilename(mfilename);
ylabel('Percentage neuronal spike coincidence with cNE events')

