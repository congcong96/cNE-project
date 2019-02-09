function [R1, R2] = ne_batch_plot_corr_measures_vs_pairwise_distance(files)

for i = 1:length(files)
    load(files{i})
    [R1{i},R2{i}] =P1figureS1(spk, spktrain, 'distance', allpwc);
end