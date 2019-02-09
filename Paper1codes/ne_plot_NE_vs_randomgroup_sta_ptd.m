function p = ne_plot_NE_vs_randomgroup_sta_ptd(files)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

stacell = cell(length(files), 1);

for i = 1:length(files)
    
    load(files{i}, 'NEgroupsta')
    stacell{i} = NEgroupsta;

end

stastruct = cell2mat(stacell');

percent = zeros(length(stastruct),1);
ptdNE = zeros(length(stastruct),1);
medsnrrandgrp = zeros(length(stastruct),1);

for i = 1:length(stastruct)
    
    staNE = stastruct(i).sta_NE;
    
    if isempty(staNE)
        continue
    end

    maxNEmat = max(staNE(:));
    minNEmat = min(staNE(:));
    ptdNE(i) = (maxNEmat - minNEmat)./stastruct(i).min_spikes;
    
    starandgrp = stastruct(i).sta_randgroups;
    starandgrp(cellfun('isempty', starandgrp)) = [];

    maxrandgrpmat = cellfun(@(x) max(x(:)), starandgrp);
    minrandgrpmat = cellfun(@(x) min(x(:)), starandgrp);
    ptdrandgrp = (maxrandgrpmat - minrandgrpmat)./stastruct(i).min_spikes;
    
    percent(i) = sum(ptdNE(i) <= [ptdNE(i);ptdrandgrp]) ./ (length(ptdrandgrp) + 1);
    
    medsnrrandgrp(i) = median(ptdrandgrp);
    
    
end

zeroidx = percent == 0;
percent(zeroidx) = [];
ptdNE(zeroidx) = [];
medsnrrandgrp(zeroidx) = [];


figure; hold on
colormap(brewermap(251,'RdYlBu'));
scatter(medsnrrandgrp, ptdNE, 10, percent)
h = refline(1,0);
h.Color = [.7 .7 .7];
h.LineStyle = '--';
tickpref;
c = colorbar;
c.TickDirection = 'out';
ylabel(c, 'Uniqueness index')

xlabel('Random group STA PTD')
ylabel('cNE STA PTD')

p(1) = signrank(ptdNE, medsnrrandgrp);
[~,p(2)] = ttest(ptdNE, medsnrrandgrp);

text(1,10,sprintf('n = %d\np < 0.001', length(ptdNE)));

print_mfilename(mfilename);

% figure;
% snrratio = snrNE./snrrandgrp;
% histogram(snrratio, 'Normalization', 'probability')
% medval = median(snrratio);
% medad = mad(snrratio, 1);
% y = ylim;
% line([medval medval], [y(1), y(2)], 'LineStyle', '--', 'Color', 'k')
% xlabel('cNE MI / neuron MI')
% ylabel('Ratio')
% tickpref;
% print_mfilename(mfilename);



