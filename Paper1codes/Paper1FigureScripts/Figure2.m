%% Plot figure A
load('site1-rn1.mat')
ne_plot_corrval_sharpness_peakdelay_graph_example(exp_site_nedata, allpwc, 15);

%% Plot figure B to G
ne_plot_PWC_statistics_and_pairwise_distance(spk, spktrain, distopt, varargin);

