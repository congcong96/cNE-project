function ne_plot_coupling_ratio_real_vs_models(realcoupratio, labels, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

assert(length(labels) == length(varargin))

L = length(labels);
modelmat = [];

for j = 1:length(varargin)
    modelmat = [modelmat varargin{j}];
end

figure;
subplot(1,2,1)
hold on
for i = 1:length(realcoupratio)
    
    plot(1:L+1, [modelmat(i,:) realcoupratio(i)], '-o','Color',[0.8,0.8,0.8],'MarkerSize',4);

end
hold off

set(gca,'xtick', 1:L+1, 'xticklabel',[labels 'real'])
xlim([0 L+2])
ylabel('Population coupling ratio')

subplot(1,2,2)
hold on
scatter(realcoupratio, varargin{3});
h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
h.LineWidth = 1.5;
hold off
xlabel('real coupling ratio')
ylabel('model coupling ratio (CR = 1)')
