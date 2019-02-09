function p = ne_plot_NE_randomgroup_sta_info(files)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

NEinfo = [];
rgrpinfo = [];
allrgrpinfo = cell(length(files),1);

for i = 1:length(files)
    
    load(files{i}, 'NEgroupinfo')
    
    allrgrpinfo{i} = {NEgroupinfo.all_randgrp_info_extrap};   
    
    NEinfo = [NEinfo NEgroupinfo.NE_info_extrap];
    rgrpinfo = [rgrpinfo NEgroupinfo.randgrp_info_extrap];


end

allrgrpinfo = [allrgrpinfo{:}];
allrgrpinfo(cellfun('isempty',allrgrpinfo)) = [];

percent = zeros(length(allrgrpinfo),1);
ctrlpercent = zeros(length(allrgrpinfo),1);


for i = 1:length(allrgrpinfo)
    
    percent(i) = sum(NEinfo(i) <= [allrgrpinfo{i};NEinfo(i)]) / (length(allrgrpinfo{i})+1);
    ctrlpercent(i) = sum(randsample(allrgrpinfo{i},1) > allrgrpinfo{i})/length(allrgrpinfo{i});
end


figure;
hold on
histogram(percent, 0:0.05:1, 'Normalization', 'probability');
x = xlim;
line([x(1) x(2)], [0.05 0.05], 'LineStyle', '--', 'Color', 'k')
% histogram(ctrlpercent, 0:0.05:1, 'Normalization', 'probability');
xlabel('p-value')
ylabel('probability')
box on
tickpref;
print_mfilename(mfilename)



figure; hold on
colormap(brewermap(251,'RdYlBu'));
scatter(rgrpinfo, NEinfo, 10, percent);

h = refline(1,0);
h.Color = [.7 .7 .7];
h.LineStyle = '--';
tickpref;


xlabel('Random group MI (bits/spike)')
ylabel('cNE MI (bits/spike)')

p(1) = signrank(NEinfo, rgrpinfo);
[~,p(2)] = ttest(NEinfo, rgrpinfo);

text(0.3,1,sprintf('n = %d\np < 0.001', length(NEinfo)));
c = colorbar;
c.TickDirection = 'out';

print_mfilename(mfilename);
