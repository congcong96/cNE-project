function ne_plot_shared_neurons_sta_corr_histogram(realcorrval, allcorrval, prctileline)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    prctileline = [];
end


if iscell(realcorrval)
    rcorrmat = cell2mat(realcorrval);
    acorrcell = cat(1,allcorrval{:});
else
    rcorrmat = realcorrval;
    acorrcell = allcorrval;
end

percent = zeros(length(rcorrmat),1);
figure;
set(gcf,'position',get(0,'screensize'))
c = 1;

for i = 1:length(rcorrmat)
    
    percent(i) = sum(rcorrmat(i) >= [acorrcell{i}; rcorrmat(i)])/ (length(acorrcell{i} + 1));
    
    subplot(4,4,c)
    hold on
    h = histogram(acorrcell{i}, 'Normalization', 'probability');
    ymax = max(h.Values);
    title(sprintf('%d',i)) % temp
%     yupp = ceil(ymax/50) * 50;
    
    line([rcorrmat(i) rcorrmat(i)], [0 ymax], 'Color','r')
%     ylim([0 yupp])
    x = xlim;
    y = ylim;
    text((x(2)-x(1))/4 + x(1), y(2)-((y(2)-y(1))/4), sprintf('p = %.3f', percent(i)));  
    
    tickpref;

    if ~isempty(prctileline)
        prc = prctile(acorrcell{i}, prctileline*100);
        line([prc prc], [0 ymax], 'Color', 'k', 'LineStyle', '--')
        legend('Simulated STA correlation', 'Real STA correlation', 'Median simulated STA correlation', 'Location', 'best')
    else
        legend('Simulated STA correlation', 'Real STRF correlation', 'Location', 'best');
    end
    if mod(c,4) == 1
        ylabel('Ratio')
    end
    if c >= 13 && c <=16
        xlabel('Correlation value')
    end
    
    hold off
    c = c+1;
    
    if mod(c,16) == 1
        figure;
        set(gcf,'position',get(0,'screensize'))
        c = 1;
    end

end

