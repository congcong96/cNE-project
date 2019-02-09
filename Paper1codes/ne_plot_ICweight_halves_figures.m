function ne_plot_ICweight_halves_figures(corrcell, shuffledcorrcell, sigvals, penetration, cNEs, prcsig)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('penetration','var')
    penetration = 11;
end

if ~exist('cNEs','var')
    cNEs = [1,1;3,3];
end

if ~exist('prcsig','var')
    prcsig = 99;
end

% plot significant correlation values histogram
figure;
edges = 0.2:0.05:1;
histogram(sigvals, edges, 'Normalization', 'probability');

text(0.3, 0.12, sprintf('n = %d', length(sigvals)));
tickpref;
print_mfilename(mfilename)

xlabel('Significant correlation values')
ylabel('Ratio')


% plot corrmat

corrmat = corrcell{penetration};
corrmat(4,:) = []; % to make figure unconfusing (quick fix)
figure; 
cmap = brewermap(1000,'reds');
colormap(cmap);
imagesc(corrmat);
c = colorbar;
c.TickDirection = 'out';
tickpref;
print_mfilename(mfilename);

xlabel('IC Bs')
ylabel('IC As')


% plot sample corrval histogram with null distribution

shuffledcorr = abs(shuffledcorrcell{penetration});
sigthresh = prctile(shuffledcorr, prcsig);

figure; hold on

histogram(shuffledcorr, 'Normalization', 'probability')

y = ylim;
line([sigthresh sigthresh], [0 y(2)], 'Color', 'b', 'LineStyle', '--')

for i = 1:size(cNEs,1)
    line([corrmat(cNEs(i,1), cNEs(i,2)) corrmat(cNEs(i,1), cNEs(i,2))], [0 y(2)])
end


legend(sprintf('Shuffled correlation values, n = %d', length(shuffledcorr)), 'Significance threshold', 'Real correlation')
tickpref;
print_mfilename(mfilename)

xlabel('Correlation value')
ylabel('Ratio')

end

