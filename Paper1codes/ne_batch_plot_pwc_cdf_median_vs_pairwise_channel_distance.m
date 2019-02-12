function [distcell, interpvalcell] = batch_ne_plot_pwc_cdf_median_vs_pairwise_channel_distance(files)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

distcell = cell(length(files),1);
interpvalcell = cell(length(files),1);
statscell = cell(length(files),1);

for i = 1:length(files)
    load(files{i})
    [statscell{i}, distcell{i}, interpvalcell{i}] = ...
        ne_plot_pwc_cdf_median_vs_pairwise_channel_distance(spk, spktrain, 'distance');
    close all
    
end


dist = cell2mat(distcell);
interpval = cell2mat(interpvalcell);

figure;
hold on
scatter(dist, interpval);

X = [ones(length(dist),1) dist];
Y = interpval;
[b, ~, ~,~, stats] = regress(Y,X);
x = min(dist)-5:max(dist)+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');


yl = ylim;
xl = xlim;

if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.8, yl(1) + (yl(2)-yl(1))*0.2, ...
        sprintf('R^{2} = %.2f\np < 0.001', stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.8, yl(1) + (yl(2)-yl(1))*0.2,...
        sprintf('R^{2} = %.2f\np = %.3f', stats(1), stats(3)))
end
tickpref;