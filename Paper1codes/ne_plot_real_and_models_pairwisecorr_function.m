function ne_plot_real_and_models_pairwisecorr_function(pwcstruct, pairs, xbounds)

realpwc = pwcstruct.real(pairs);
shuffpwc = pwcstruct.shuffled(pairs);
PRpwc = pwcstruct.PR(pairs);
DGpwc = pwcstruct.DG(pairs);

for j = 1:length(realpwc)
    
    figure;
    
    subplot(2,2,1) 
    delay = realpwc(j).delay;
    r12 = realpwc(j).r12;
    hold on;
    bar(delay, r12, 'k');
    tickpref;
    title(sprintf('Neuron %d and %d',...
        realpwc(j).pairs(1),realpwc(j).pairs(2)));
    if ~isempty(xbounds)
        xlim([-xbounds xbounds])
    end
    y = ylim;
    
    subplot(2,2,2) 
    delay = shuffpwc(j).delay;
    r12 = shuffpwc(j).r12;
    hold on;
    bar(delay, r12, 'k');
    tickpref;
    title(sprintf('Neuron %d and %d',...
        shuffpwc(j).pairs(1),shuffpwc(j).pairs(2)));
    if ~isempty(xbounds)
        xlim([-xbounds xbounds])
    end
    ylim([y(1) y(2)])
    
    subplot(2,2,3) 
    delay = PRpwc(j).delay;
    r12 = PRpwc(j).r12;
    hold on;
    bar(delay, r12, 'k');
    tickpref;
    title(sprintf('Neuron %d and %d',...
        PRpwc(j).pairs(1),PRpwc(j).pairs(2)));
    if ~isempty(xbounds)
        xlim([-xbounds xbounds])
    end
    ylim([y(1) y(2)])
    
    subplot(2,2,4) 
    delay = DGpwc(j).delay;
    r12 = DGpwc(j).r12;
    hold on;
    bar(delay, r12, 'k');
    tickpref;
    title(sprintf('Neuron %d and %d',...
        DGpwc(j).pairs(1),DGpwc(j).pairs(2)));
    if ~isempty(xbounds)
        xlim([-xbounds xbounds])
    end
    ylim([y(1) y(2)])
    print_mfilename(mfilename);

end



