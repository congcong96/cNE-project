function thresh = ne_plot_spon_evoked_corrval_histogram(corrvalstruct, varargin)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

assert(length(corrvalstruct) == 1, 'Choose only one penetration to run!')

nulldist = corrvalstruct.shuffledcorrvals;
maxnulldist = max(abs(nulldist));
edges = 0:0.02:maxnulldist;


sponcorr = corrvalstruct.sponcorrvals;
evokedcorr = corrvalstruct.evokedcorrvals;

figure;
c = 1;

for i = 1:length(sponcorr)
    subplot(3,2,c)
    hold on
    histogram(abs(nulldist), edges, 'Normalization', 'probability')
    xlabel('Correlation value')
    ylabel('Ratio')
    y = ylim;
    sponval = abs(sponcorr(i));
    line([sponval sponval], [0 y(2)], 'Color','r','LineStyle','--')
    if corrvalstruct.sponpvals(i) < 0.0005
        text(0.5, 0.1, sprintf('p = %.3e', corrvalstruct.sponpvals(i)));
    else
        text(0.5, 0.1, sprintf('p = %.3f', corrvalstruct.sponpvals(i)));
    end
    legend(sprintf('Shuffled correlation values, n = %d', length(nulldist)), 'Real correlation value')
    tickpref;
    
    if c == 6 && i ~= length(sponcorr)
        print_mfilename(mfilename)
        c = 1;
        figure;
    else
        c = c+1;
    end
    
end

print_mfilename(mfilename)

figure;
c = 1;

for i = 1:length(evokedcorr)
    subplot(3,2,c)
    hold on
    histogram(abs(nulldist), edges, 'Normalization', 'probability')
    xlabel('Correlation value')
    ylabel('Ratio')
    y = ylim;
    evokedval = abs(evokedcorr(i));
    line([evokedval evokedval], [0 y(2)], 'Color','r','LineStyle','--')
    if corrvalstruct.evokedpvals(i) < 0.0005
        text(0.5, 0.1, sprintf('p = %.3e', corrvalstruct.evokedpvals(i)));
    else
        text(0.5, 0.1, sprintf('p = %.3f', corrvalstruct.evokedpvals(i)));
    end
    legend(sprintf('Shuffled correlation values, n = %d', length(nulldist)), 'Real correlation value')
    tickpref;
    
    if c == 6 && i ~= length(evokedcorr)
        print_mfilename(mfilename)
        c = 1;
        figure;
    else
        c = c+1;
    end
    
end

print_mfilename(mfilename)



if ~isempty(varargin)
    
    figure;
    hold on
    histogram(abs(nulldist),edges,'Normalization','probability')
    thresh = prctile(abs(nulldist),99);
    y = ylim;
    line([thresh thresh], [0 y(2)], 'Color', 'b', 'LineStyle', '--')
    
    for i = 1:length(varargin{1})
        line([abs(sponcorr(varargin{1}(i))) abs(sponcorr(varargin{1}(i)))], [0 y(2)])
    end
       
    legend(sprintf('Shuffled correlation values, n = %d', length(nulldist)), 'Significance threshold', 'Real correlation')
    tickpref;
    
end
    

