function ne_scroll_through_NE_activity_graph(exp_site_nedata, startx)

% startx in seconds
% Written 11/3/16 by JS.

if ~exist('startx', 'var') || isempty(startx)
    startx = 0;
end

% plot all NEactivity in one graph using distinguishable colors
% NEactivity = exp_site_nedata.nedata.Activities;
NEactivity = exp_site_nedata.nedata.total_activities;
figure;
hold on
linecolors = distinguishable_colors(size(NEactivity,1));
leg = cell(size(NEactivity,1),1);
for i = 1:size(NEactivity,1)
    plot(NEactivity(i,:), 'Color', linecolors(i,:));
    leg{i} = num2str(i);
end

% set ymax to 95th percentile of all activity
ymax = prctile(NEactivity(:), 99.5);

xlabel('Time (s)', 'FontSize', 14);
tickpref;    
legend(leg,'Location','Best');

% to label x axis correctly
df = exp_site_nedata.df;
binsize = df * 0.5; % binsize in ms
snippetlength = 200*binsize; %length of snippet in ms
xtick = 0:50:size(NEactivity,2);
xticklabel = 0:snippetlength/4/1000:binsize/1000*size(NEactivity,2);
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
set(gcf,'position', [432 334 1121 620]);
box on

hold off

%convert startx to sample number
startx = startx*1000/binsize;

while startx <= size(NEactivity,2) - 200
    xlim([startx startx + 200])
%     ylim([-20 ymax])
    
    startx = startx + 100;
    pause;
end