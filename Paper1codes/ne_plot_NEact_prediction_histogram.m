function ne_plot_NEact_prediction_histogram(corractmat, sortopt)

if nargin == 1 
    sortopt = 0;
end

if sortopt == 1
    
    offdiags = cell2mat(cellfun(@(x) x(~eye(size(x,1))), corractmat, 'UniformOutput', 0));
    minoffdiag = floor(min(offdiags)*10)/10; % get nearest tenth
    diags = cell2mat(cellfun(@diag, corractmat, 'UniformOutput', 0));

    figure;
    hold on
    histogram(offdiags, minoffdiag:0.05:1, 'FaceColor', 'b')% 'Normalization', 'probability');
    histogram(diags, minoffdiag:0.05:1, 'FaceColor', 'r')% 'Normalization', 'probability');
    ylabel('Count')
    xlabel('NE activity correlation')
    tickpref;
    print_mfilename(mfilename)

    % plot 95th percentile of offdiags
    prc = prctile(offdiags, 95);
    y = ylim;
    line([prc prc], [y(1) y(2)], 'Color', 'b', 'LineStyle', '--');

    legend(sprintf('non-matching, n = %d',length(offdiags)),...
        sprintf('matching, n = %d', length(diags)),...
        'non-matching 95th percentile','Location','Best');
    
    
else
    
    allcorrvals = cell2mat(cellfun(@(x) x(:), corractmat, 'UniformOutput', 0));
    expected_matches = sum(cellfun(@(x) min(size(x)), corractmat));
    totalcorrvals = length(allcorrvals);
    ratio = totalcorrvals/expected_matches;
    mincorrval = floor(min(allcorrvals)*10)/10;
    
    figure; 
    histogram(allcorrvals, mincorrval:0.05:1, 'Normalization', 'probability');
    xlabel('NE activity correlation')
    
    ylab = get(gca, 'YTickLabel');
    ytick = get(gca, 'YTick');
    normalizedright = cell2mat(cellfun(@(x) str2double(x).*ratio, ylab, 'UniformOutput', 0));    
    rightsteps = 0:0.1:max(normalizedright);
    rightlabels = interp1(normalizedright, ytick, rightsteps);
    ylabel('ratio of all correlation values')        
    
    yyaxis right
    histogram(allcorrvals, mincorrval:0.05:1, 'Normalization', 'probability', 'FaceColor', 'b');
    set(gca, 'YTick', rightlabels, 'YTickLabel', rightsteps);
    ylabel('ratio of expected matches (number of cNEs)')
    
    tickpref;
    print_mfilename(mfilename)
    
end


