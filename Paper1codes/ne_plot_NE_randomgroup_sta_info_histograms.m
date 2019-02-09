function percent = ne_plot_NE_randomgroup_sta_info_histograms(NEgroupinfo)

c = 1;
figure;

percent = zeros(length(NEgroupinfo),1);

for i = 1:length(NEgroupinfo)
    
    NEinfo = NEgroupinfo(i).NE_info_extrap;
    randgrpinfo = NEgroupinfo(i).all_randgrp_info_extrap;
    medpNEinfo = median(randgrpinfo);
    percent(i) = 1 - (sum(NEinfo > [randgrpinfo; NEinfo])/ (length(randgrpinfo)+1) );
    
    if isempty(NEinfo)
        continue
    end
    
    subplot(3,2,c); hold on;
    histogram(randgrpinfo, 'Normalization', 'probability')
    y = ylim;
    line([medpNEinfo medpNEinfo], [y(1) y(2)], 'Color', 'b', 'LineStyle', '--')
    line([NEinfo NEinfo], [y(1) y(2)], 'Color', 'r', 'LineStyle', '--')
    legend('random group MI (n = 500)', 'median random group MI', 'cNE MI', 'Location', 'best')
    x = xlim;
    text(x(2) - (x(2) - x(1))/4, y(2) - (y(2) - y(1))/2, sprintf('p = %.3f', percent(i)));
    tickpref;
    
    c = c+1;
    
    if c == 7 && c <= length(NEgroupinfo)
        print_mfilename(mfilename);

        figure;
        c = 1;
    end

end

print_mfilename(mfilename);