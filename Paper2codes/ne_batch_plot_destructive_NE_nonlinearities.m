function ne_batch_plot_destructive_NE_nonlinearities(fiofit)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(fiofit)
    
    figure;
    
    for j = 1:length(fiofit{i})
        
        subplot(2,3,j)
        hold on
        plot(fiofit{i}(j).bins, fiofit{i}(j).pspkx, 'ko', 'markerfacecolor', 'k');
        plot(fiofit{i}(j).xFit, fiofit{i}(j).pspkxFit, 'r-');
        xmax = max([max(fiofit{i}(j).bins) abs(min(fiofit{i}(j).bins))]);
        xmax = xmax + 2*xmax*0.05;
        xlim([-xmax xmax]);
        plot([-xmax xmax], [fiofit{i}(j).pspk fiofit{i}(j).pspk], 'k--');
        ylimit = get(gca,'ylim');
        set(gca,'ylim', [-0.1*max(fiofit{i}(j).pspkx) max(ylimit)]);
        set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
        legend('Data', 'HVV Fit', 'location', 'northwest');
        title(sprintf('NMSE = %.3f, R2 = %.3f', fiofit{i}(j).nmse, fiofit{i}(j).r2));
        
    end
    


end

