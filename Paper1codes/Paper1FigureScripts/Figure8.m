%% Figure 8D

% Load example penetration
evoked = {'170823_212426-site3-1402um-30db-rn1-30min-H31x64-fs20000-spk-filt-ne-20dft.mat'};
spon = {'170823_212426-site3-1402um-30db-spon-30min-H31x64-fs20000-spk-filt-ne-20dft.mat'};

% Get correlation matrix and plot it
[corrmat, orderidx] = ne_batch_plot_arranged_ICweight_corrmat(spon, evoked, 'max_row');


%% Figure 8B

% Get example ICs for stem plots
spondata = load(spon{1});
evokeddata = load(evoked{1});
sponIC = spondata.exp_site_nedata.nedata.Patterns;
evokedIC = evokeddata.exp_site_nedata.nedata.Patterns;

% Define stem plot colors
color = eight_color_blind_palette('blue', 'vermillion');

ne_plot_ICweights_of_halves_examples(sponIC, evokedIC, orderidx{1}, color);

%% Figure 8C

% Get pre-calculated correlation values for shuffled data
load('corrvals_112618.mat')
idx = 7;
histcolors = fifteen_color_blind_palette(10, 14, 2);
thresh = prctile(abs(corrvalstruct(idx).shuffledcorrvals), 99);

figure; hold on

% Plot histogram of shuffled correlation values
histogram(abs(corrvalstruct(idx).shuffledcorrvals), 0:0.02:0.8, 'Normalization', 'probability', 'FaceColor', histcolors(1,:));
y = ylim;

% Plot significance threshold
line([thresh thresh], [y(1) y(2)], 'Color', histcolors(1,:), 'LineStyle', '--');

% Plot real correlation values of examples from Figure 8B
line([abs(corrvalstruct(idx).evokedcorrvals(6)) abs(corrvalstruct(idx).evokedcorrvals(6))],...
    [y(1) y(2)], 'Color', histcolors(2,:), 'LineStyle', '-');
line([abs(corrvalstruct(idx).evokedcorrvals(10)) abs(corrvalstruct(idx).evokedcorrvals(10))],...
    [y(1) y(2)], 'Color', histcolors(3,:), 'LineStyle', '-');

% Make the plot pretty 
xlim([0 0.9])
legend('Shuffled correlation, n = 100000', 'Significance threshold (p = 0.01)', 's-cNE 9 vs e-cNE 6', 's-cNE 5 vs e-cNE 10', 'Location', 'best');
tickpref;
ylabel('Ratio')
xlabel('|Correlation|')

%% Figure 8E

[perc, sigcorrvals] = ne_calc_spon_evoked_percentage_matches(corrvalstruct);
