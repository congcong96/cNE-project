% load('150722_225641-site5-812um-20db-rn16-fs20000-A-spk-strfcmb-ne-20dft.mat', 'exp_site_nedata')
% load('NEacthalves.mat', 'shuffledactcorrcell')
% [actcorr, trainingact, testact] = ne_predict_subset_of_NE_activity_with_ICweights(exp_site_nedata, 'plotopt', 1, 'plottime', 0);
% 
% figure;
% cmap = brewermap(1000,'reds');
% imagesc(actcorr);
% colormap(cmap);
% c = colorbar;
% c.TickDirection = 'out';
% tickpref;
% print_mfilename(mfilename);
% 
% figure;
% shuffcorrvec = cell2mat(shuffledactcorrcell);
% histogram(shuffcorrvec, 'Normalization', 'probability');
% sigthresh = prctile(shuffcorrvec, 99.9);
% 
% y = ylim;
% line([sigthresh sigthresh], [0 y(2)], 'Color', 'b', 'LineStyle', '--')
% line([0.85 0.85], [0 y(2)], 'Color', 'g', 'LineStyle', '-')
% line([0.99 0.99], [0 y(2)], 'Color', 'k', 'LineStyle', '-')
% print_mfilename(mfilename)

function perc = plot_figure4s3(actcorrcell, prcthresh)

if ~exist('prcthresh','var')
    prcthresh = 99;
end

maxactrow = cellfun(@(x) max(x, [], 2), actcorrcell, 'UniformOutput', 0);
nonmatches = cellfun(@(x,y) setdiff(x,y), actcorrcell, maxactrow, 'UniformOutput', 0);

% plot matched vs non-matched histogram
figure;

% yyaxis left
nonmatches = cell2mat(nonmatches);
histogram(nonmatches, -0.2:0.05:1, 'Normalization', 'probability')
ylabel(sprintf('Ratio (non-matched), n = %d', length(nonmatches)))

% yyaxis right
% maxactrow = cell2mat(maxactrow);
% histogram(maxactrow, -0.2:0.05:1, 'Normalization', 'probability')
% ylabel(sprintf('Ratio (matched), n = %d', length(maxactrow)))

xlabel('Correlation')

tickpref;
sigthresh = prctile(nonmatches, prcthresh);
y = ylim;
line([sigthresh sigthresh], [0 y(2)], 'Color', 'b', 'LineStyle', '--')
samplelow = actcorrcell{12}(4,6);
samplehigh = actcorrcell{12}(6,5);
line([samplelow samplelow], [0 y(2)], 'Color', 'g', 'LineStyle', '-')
line([samplehigh samplehigh], [0 y(2)], 'Color', 'k', 'LineStyle', '-')
legend('Non-matched correlations', 'Significance threshold', 'Sample 1', 'Sample 2')
print_mfilename(mfilename);

% plot corrmat
corrmat = actcorrcell{12};
corrmat(:,4) = []; % to make figure unconfusing (quick fix)
figure; 
cmap = brewermap(1000,'reds');
colormap(cmap);
imagesc(corrmat);
c = colorbar;
c.TickDirection = 'out';
tickpref;
print_mfilename(mfilename);

xlabel('Test cNE activity')
ylabel('Training cNE activity')

% Calculate percentage significant
maxactrow = cell2mat(maxactrow);
perc = sum(maxactrow > sigthresh) / length(maxactrow);

figure;
sigactrow = maxactrow(maxactrow >= sigthresh);
histogram(sigactrow, 'Normalization', 'probability')
ylabel(sprintf('Ratio (matched), n = %d', length(sigactrow)))
y = ylim;
% line([sigthresh sigthresh], [0 y(2)], 'Color', 'b', 'LineStyle', '--')
print_mfilename(mfilename);
tickpref;
xlabel('Correlation')





end
