function ne_plot_NE_vs_randomgroup_sta_snr_histogram_example(NEgroupsta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

percent = zeros(length(NEgroupsta),1);
ptdNE = zeros(length(NEgroupsta),1);
medptdrandgrp = zeros(length(NEgroupsta),1);

c = 1;
figure;

for i = 1:length(NEgroupsta)
    
    staNE = NEgroupsta(i).sta_NE;
    
    if isempty(staNE)
        continue
    end

    maxNEmat = max(staNE(:));
    minNEmat = min(staNE(:));
    ptdNE(i) = (maxNEmat - minNEmat)./NEgroupsta(i).min_spikes;
    
    starandgrp = NEgroupsta(i).sta_randgroups;
    starandgrp(cellfun('isempty', starandgrp)) = [];

    maxrandgrpmat = cellfun(@(x) max(x(:)), starandgrp);
    minrandgrpmat = cellfun(@(x) min(x(:)), starandgrp);
    ptdrandgrp = (maxrandgrpmat - minrandgrpmat)./NEgroupsta(i).min_spikes;
    
    percent(i) = sum(ptdNE(i) <= [ptdNE(i);ptdrandgrp]) ./ (length(ptdrandgrp) + 1);
    
    medptdrandgrp(i) = median(ptdrandgrp);
    
    subplot(3,2,c)
    histogram(ptdrandgrp, 'Normalization', 'probability')
    y = ylim;
    line([medptdrandgrp(i) medptdrandgrp(i)], [y(1) y(2)], 'Color', 'b', 'LineStyle', '--')
    line([ptdNE(i) ptdNE(i)], [y(1) y(2)], 'Color', 'r', 'LineStyle', '--')
    legend('random group MI (n = 500)', 'median random group MI', 'cNE MI', 'Location', 'best')
    x = xlim;
    text(x(2) - (x(2) - x(1))/4, y(2) - (y(2) - y(1))/2, sprintf('p = %.3f', percent(i)));
    tickpref;
    
    if c == 6 && i ~= length(NEgroupsta)
        print_mfilename(mfilename);

        figure;
        c = 1;
    else
        c = c+1;
    end    
    
end
print_mfilename(mfilename);

