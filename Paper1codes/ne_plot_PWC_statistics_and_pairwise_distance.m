function [R1, R2] = ne_plot_PWC_statistics_and_pairwise_distance(spk, spktrain, distopt, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fsd = 2000;
stimdur = 1/fsd*1000 * size(spktrain,2);

if nargin == 4
    pwc = varargin{1};
elseif nargin < 4
    pwc = pairwisecorr_function(spktrain, stimdur, fsd);
else
    error('Too many inputs')
end

[~, ~, interpvals, ~, pwcidx] = calc_pwc_sharpness_cdf(pwc, 40, 0.5, 0);
interpvals(isnan(interpvals)) = 0;

newpwc = pwc(pwcidx);

peakdelay = abs([newpwc.peakdelay]');
pairs = cell2mat({newpwc.pairs}');

dsspktrain = downsample_spiketrain(spktrain, 20);
corrmat = corr(dsspktrain');

corrvec = arrayfun(@(x,y) corrmat(x,y),pairs(:,1),pairs(:,2));


% get pairs and position
pairs = cell2mat({pwc(pwcidx).pairs}');
pos = cell2mat({spk.position}');

% calculate distance
switch distopt
    case 'distance'
        dist = ne_calc_interneuronal_distances(pairs, pos, 'distance');
    case 'channel'
        dist = ne_calc_interneuronal_distances(pairs, pos, 'channel');
end


% graphtitle = sprintf('%s-site%d-%s',exp_site_nedata.exp, exp_site_nedata.site, exp_site_nedata.stim);
% %correct underscores
% graphtitle = regexprep(graphtitle, '_', '\\_');

figure;
set(gcf, 'Position', [200 50 1500 900])

subplot(3,2,1)
% posidx = corrvec >= 0;
% [f, gof{1}] = fit(log(corrvec(posidx)),interpvals(posidx),'exp1');
[f, gof{1}] = fit(corrvec,interpvals,'exp1');
% plot(f,log(corrvec(posidx)),interpvals(posidx))
plot(f,corrvec,interpvals)

xlabel('spike train Pearson''s correlation')
ylabel('PWC sharpness')
% title(graphtitle)

x = xlim;
y = ylim;
text((x(2) - x(1))/4*3+x(1), (y(2) - y(1)) /5*3+y(1), sprintf('R^{2} = %.3f', gof{1}.adjrsquare));  
tickpref;

lgd = legend;
lgd.Location = 'best';
R1(1) = gof{1}.adjrsquare;

subplot(3,2,3)
% [f, gof{2}] = fit(log(corrvec(posidx)),peakdelay(posidx),'exp1');
[f, gof{2}] = fit(corrvec,peakdelay,'exp1');
% plot(f,log(corrvec(posidx)),peakdelay(posidx))
plot(f,corrvec,peakdelay)

xlabel('spike train Pearson''s correlation')
ylabel('PWC peak delay')
% title(graphtitle)

x = xlim;
y = ylim;
text((x(2) - x(1))/4*3+x(1), (y(2) - y(1)) /5*3+y(1), sprintf('R^{2} = %.3f', gof{2}.adjrsquare));  
tickpref;

lgd = legend;
lgd.Location = 'best';
R1(2) = gof{2}.adjrsquare;


subplot(3,2,5)
[f, gof{3}] = fit(interpvals, peakdelay, 'exp1');
plot(f, interpvals, peakdelay)

xlabel('PWC sharpness')
ylabel('PWC peak delay')
% title(graphtitle)

x = xlim;
y = ylim;
text((x(2) - x(1))/6+x(1), (y(2) - y(1)) /5*3+y(1), sprintf('R^{2} = %.3f', gof{3}.adjrsquare));  
tickpref;

lgd = legend;
lgd.Location = 'best';
R1(3) = gof{3}.adjrsquare;


%remove dist == 0
zeroidx = dist == 0;
dist(zeroidx) = [];
interpvals(zeroidx) = [];
peakdelay(zeroidx) = [];
corrvec(zeroidx) = [];


subplot(3,2,2)
hold on
scatter(dist, interpvals, 30, '.');

X = [ones(length(dist),1) dist];
Y = interpvals;
[b, ~, ~,~, stats] = regress(Y,X);
x = min(dist)-5:max(dist)+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');


yl = ylim;
xl = xlim;

if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.6, yl(1) + (yl(2)-yl(1))*0.2, ...
        sprintf('R^{2} = %.2f\np < 0.001', stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.6, yl(1) + (yl(2)-yl(1))*0.2,...
        sprintf('R^{2} = %.2f\np = %.3f', stats(1), stats(3)))
end
tickpref;
xlabel('distance (\mum)')
ylabel('sharpness')
R2(1) = stats(1);

subplot(3,2,4)
hold on
scatter(dist, peakdelay, 30,'.');

X = [ones(length(dist),1) dist];
Y = peakdelay;
[b, ~, ~,~, stats] = regress(Y,X);
x = min(dist)-5:max(dist)+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');


yl = ylim;
xl = xlim;

if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.6, yl(1) + (yl(2)-yl(1))*0.2, ...
        sprintf('R^{2} = %.2f\np < 0.001', stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.6, yl(1) + (yl(2)-yl(1))*0.2,...
        sprintf('R^{2} = %.2f\np = %.3f', stats(1), stats(3)))
end
tickpref;
xlabel('distance (\mum)')
ylabel('peak delay (ms)')
R2(2) = stats(1);


subplot(3,2,6)
hold on

scatter(dist, corrvec, 30, '.');

X = [ones(length(dist),1) dist];
Y = corrvec;
[b, ~, ~,~, stats] = regress(Y,X);
x = min(dist)-5:max(dist)+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');


yl = ylim;
xl = xlim;

if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.6, yl(1) + (yl(2)-yl(1))*0.2, ...
        sprintf('R^{2} = %.2f\np < 0.001', stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.6, yl(1) + (yl(2)-yl(1))*0.2,...
        sprintf('R^{2} = %.2f\np = %.3f', stats(1), stats(3)))
end
tickpref;
xlabel('distance (\mum)')
ylabel('correlation value')
R2(3) = stats(1);


end

