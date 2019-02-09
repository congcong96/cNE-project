function [NEhw, nNEhw, varargout] = ne_calc_NE_vs_nonNE_strf_halfwidth(exp_site_nedata, strf, trigger, plotopt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    plotopt = 0;
end


%get cNE pairs
pwcNE = exp_site_nedata.nedata.pwc;
pwcNE = cell2mat(pwcNE);
NEpairs = cell2mat({pwcNE.pairs}');
[~, uniqueidx] = unique(NEpairs, 'rows', 'stable');
NEpairs = NEpairs(uniqueidx,:);
NEhw = zeros(size(NEpairs,1),1);

%get non-cNE pairs
allpairs = nchoosek(1:length(strf), 2);
[~, nonNEidx] = setdiff(allpairs, NEpairs, 'rows');
nonNEpairs = allpairs(nonNEidx,:);
nNEhw = zeros(size(nonNEpairs,1),1);


for i = 1:size(NEpairs,1)
    neuron1 = NEpairs(i,1);
    neuron2 = NEpairs(i,2);
    [c12, lags] = calc_STRF_xcorr(strf, trigger, neuron1, neuron2);
    NEhw(i) = find_strf_corr_halfwidth(lags,c12);
end

for i = 1:size(nonNEpairs,1)
    neuron1 = nonNEpairs(i,1);
    neuron2 = nonNEpairs(i,2);
    [c12, lags] = calc_STRF_xcorr(strf, trigger, neuron1, neuron2);
    nNEhw(i) = find_strf_corr_halfwidth(lags,c12);
end


if plotopt == 1
  
%     figure;
    boxinput = [NEhw;nNEhw];
    group = [zeros(length(NEhw),1);ones(length(nNEhw),1)];

    plotSpread(NEhw, 'distributionColors', [1 .6 .6], 'xValues', 1, 'spreadWidth', 0.5)
    plotSpread(nNEhw, 'distributionColors', [.6 .6 1], 'xValues', 2, 'spreadWidth', 0.5)

    boxplot(boxinput, group, 'notch','on','labels',{'NE', 'nonNE'},...
        'colors','rb', 'symbol', 'b','outliersize',3, 'widths', 0.2);

%     ylabel('STRF cross-correlation half-width (ms)')
    
    p = ranksum(NEhw, nNEhw);
    y = ylim;
    text(1.3, y(1) + ((y(2)-y(1))/2), sprintf('p = %.3f',p));
    tickpref;

%     print_mfilename(mfilename);
    
    varargout{1} = p;
    

end

end

