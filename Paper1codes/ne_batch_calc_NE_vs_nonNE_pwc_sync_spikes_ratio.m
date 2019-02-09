function ne_batch_calc_NE_vs_nonNE_pwc_sync_spikes_ratio(NEfiles, plotopt)

NEratio = cell(length(NEfiles),1);
nonNEratio = cell(length(NEfiles),1);

for i = 1:length(NEfiles)
    load(NEfiles{i})
    basefile = regexp(NEfiles{i}, '^\S+(?=(-ne))','match','once');
    spkfile = [basefile '.mat'];
    load(spkfile, 'allpwc')
    [NEratio{i}, nonNEratio{i}] = ne_calc_NE_vs_nonNE_pwc_sync_spikes_ratio(exp_site_nedata, allpwc);
end

if plotopt == 1
    NEmat = cell2mat(NEratio');
    nonNEmat = cell2mat(nonNEratio');
    NEedges = 0:0.05:max(NEmat);
    nonNEedges = 0:0.05:max(nonNEmat);
    
    figure;
    histogram(NEmat, NEedges, 'Normalization', 'probability');
    xlabel('Ratio of synchronized spikes')
    ylabel('Probability')
    title('NE pairs')
    tickpref;
    x = xlim;
    y = ylim;
    text(x(2) - x(2)/4, y(2) - y(2)/4, sprintf('n = %d', length(NEmat)));
    print_mfilename(mfilename);
    
    figure;
    histogram(nonNEmat, nonNEedges, 'Normalization', 'probability');
    xlabel('Ratio of synchronized spikes')
    ylabel('Probability')
    title('non-NE pairs')
    tickpref;
    x = xlim;
    y = ylim;
    text(x(2) - x(2)/4, y(2) - y(2)/4, sprintf('n = %d', length(nonNEmat)));
    print_mfilename(mfilename);
    
    figure; hold on
    plotSpread(NEmat', 'distributionColors', [.8 .8 .8], 'xValues',...
        1, 'spreadWidth', .5);
    plotSpread(nonNEmat', 'distributionColors', [.8 .8 .8], 'xValues',...
        2, 'spreadWidth', .5);
    
    boxinput = [NEmat, nonNEmat];
    bsvec = [ones(1,length(NEmat)) 2 * ones(1,length(nonNEmat))];
    
    boxplot(boxinput, bsvec, 'notch','on', 'colors','b', 'symbol',...
        'b','outliersize',3, 'widths', .5, 'labels', {'NE ratio','non-NE ratio'})  
    
    ylabel('Ratio of synchronized spikes')
    tickpref;
    [~, p] = ttest2(NEmat, nonNEmat);
    x = xlim;
    y = ylim;
    text(x(1) + (x(2) - x(1))/2, y(2)- (y(2)-y(1))/4, sprintf('p = %.3e', p));
    
    print_mfilename(mfilename);

   
end


    