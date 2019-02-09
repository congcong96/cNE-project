function ne_plot_ICweights_of_halves_examples(IC1, IC2, idx, color)

if isempty(idx)
    corrmat = abs(corr(IC1, IC2));    
    maxcol = max(corrmat, [], 1);
    sigcols = (maxcol > 0.30);
    threshpatcorr = corrmat(:, sigcols);
    [~, maxidx] = max(threshpatcorr, [], 1); 
    sortedIC2 = IC2(:,sigcols);
    sortedIC1 = IC1(:,maxidx);
else
    sortedIC1 = IC1(:,idx);
    sortedIC2 = IC2;
end

diagonal = diag(corr(sortedIC1, sortedIC2));

for i = 1:length(diagonal)
    if diagonal(i) < 0
        sortedIC2(:,i) = -sortedIC2(:,i);
    end
    
    ymax = max(abs([sortedIC1(:,i); sortedIC2(:,i)])) + 0.05;
    
    figure;
    
    if exist('color','var')
        set(gcf,'defaultAxesColorOrder', color);
    end
    
    hold on
    yyaxis left
    stem(sortedIC1(:,i))
    ylim([-ymax ymax])
    ylabel('Independent component weight', 'Color', 'k')

    
    yyaxis right
    axis ij
    stem(sortedIC2(:,i))
    ylim([-ymax ymax])
    
    legend('IC A', 'IC B', 'Location', 'Best')
    xlabel('Neuron number')
%     ylabel('Independent component weight')
    
    xlim ([0 size(sortedIC1,1) + 1]);
    
    x = xlim;
    y = ylim;
    
    text(x(2)/5, y(2) - (y(2)-y(1))/5, sprintf('|correlation| = %.2f', abs(diagonal(i))))
    
    tickpref;
    print_mfilename(mfilename)
    
end
