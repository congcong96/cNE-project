function ne_plot_strfcorr_vs_pwc_halfwidth_comparison(hw)

% plot hw scatter plot
figure;
hold on
scatter([hw.nonNE.strfhw], [hw.nonNE.cchw], 10, 'b', '.');
scatter([hw.NE.strfhw], [hw.NE.cchw],10, 'r', '.');
h = refline(1,0);
h.Color = 'k';
h.LineStyle = '--';

% plot cchw / strfhw
figure;
NEratio = [hw.NE.cchw]./[hw.NE.strfhw];
nonNEratio = [hw.nonNE.cchw]./[hw.nonNE.strfhw];

NEratio(isnan(NEratio)) = [];
nonNEratio(isnan(nonNEratio)) = [];

plotSpread(NEratio', 'distributionColors', [1 .6 .6], 'xValues',...
    1, 'spreadWidth', .5);
plotSpread(nonNEratio', 'distributionColors', [.6 .6 1], 'xValues',...
    2, 'spreadWidth', 0.5)

% mnNEratio = mean(NEratio);
% mnnonNEratio = mean(nonNEratio);
% semNEratio = std(NEratio) / sqrt(length(NEratio));
% semnonNEratio = std(nonNEratio) / sqrt(length(nonNEratio));
% errorbar([1 2], [mnNEratio mnnonNEratio], [semNEratio semnonNEratio]);
boxinput = [NEratio nonNEratio];
group = [zeros(length(NEratio),1);ones(length(nonNEratio),1)];

boxplot(boxinput, group, 'notch','on','labels',{'NE', 'nonNE'},...
    'colors','rb', 'symbol', 'b','outliersize',3, 'widths', 0.2);



