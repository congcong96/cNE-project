function ne_plot_sta_corrvals(labels, varargin)

boxinput = [];
boxlabel = [];
figure;
hold on

for i = 1:length(labels)
    plotSpread(varargin{i}, 'distributionColors', [.8 .8 .8], 'xValues',...
        i, 'spreadWidth', 0.5)
    
    boxinput = [boxinput; varargin{i}];
    boxlabel = [boxlabel; i * ones(length(varargin{i}))];
end

boxplot(boxinput, boxlabel, 'notch','on','labels', labels, 'colors','b',...
    'symbol', 'b','outliersize',3, 'widths', .5)

ylabel('STA correlation values')
tickpref;
print_mfilename(mfilename);