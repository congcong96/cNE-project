function varargout = ne_plot_pwc_sharpness_peakdelay_vs_pairwise_channel_distance(spk, spktrain, distopt, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%initialize values
fsd = 2000;
stimdur = 1/fsd*1000 * size(spktrain,2);

%calculate PWCs if not already calculated
if nargin == 4
    pwc = varargin{1};
elseif nargin < 4
    pwc = pairwisecorr_function(spktrain, stimdur, fsd);
else
    error('Too many inputs')
end

%calculate sharpness
[~, ~, interpval, pwcidx] = calc_pwc_sharpness_cdf(pwc, 40, 0.5, 0);
interpval(isnan(interpval)) = 0;

% get pairs and position
pairs = cell2mat({pwc(pwcidx).pairs}');
pos = cell2mat({spk.position}');

% calculate distance
switch distopt
    case 'distance'
        dist = ca_calc_interneuronal_distances(pairs, pos, 'distance');
    case 'channel'
        dist = ca_calc_interneuronal_distances(pairs, pos, 'channel');
end

% get peakdelay
peakdelay = [pwc(pwcidx).peakdelay];

% calculate Pearson's corr at 10ms time bins
dsspktrain = downsample_spiketrain(spktrain, 20);
corrmat = corr(dsspktrain');
corrvec = arrayfun(@(x,y) corrmat(x,y), pairs(:,1),pairs(:,2));


%remove dist == 0
zeroidx = dist == 0;
dist(zeroidx) = [];
interpval(zeroidx) = [];
peakdelay(zeroidx) = [];
corrvec(zeroidx) = [];


figure;
subplot(2,2,1)
hold on
scatter(dist, interpval);

X = [ones(length(dist),1) dist];
Y = interpval;
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

subplot(2,2,2)
hold on
scatter(dist, peakdelay);

X = [ones(length(dist),1) dist];
Y = peakdelay';
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

subplot(2,2,3)
hold on

scatter(dist, corrvec);

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


temp = {stats dist interpval peakdelay pwc};


for i = 1:nargout
    varargout{i} = temp{i};
end

