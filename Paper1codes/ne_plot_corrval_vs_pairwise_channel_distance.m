function varargout = ne_plot_corrval_vs_pairwise_channel_distance(exp_site_nedata, distopt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% 
% fsd = 2000;
% stimdur = 1/fsd*1000 * size(spktrain,2);
% 
% if nargin == 4
%     pwc = varargin{1};
% elseif nargin < 4
%     pwc = pairwisecorr_function(spktrain, stimdur, fsd);
% else
%     error('Too many inputs')
% end
% 
% [~, ~, interpval, pwcidx] = calc_pwc_sharpness_cdf(pwc, 40, 0.5, 0);
% 
% pairs = cell2mat({pwc(pwcidx).pairs}');
% 
% pos = cell2mat({spk.position}');

pos = cell2mat(cellfun(@str2num, exp_site_nedata.nedata.position, 'UniformOutput', 0));
pairs = nchoosek(1:size(pos,1), 2);

corrmat = corr(exp_site_nedata.nedata.spktrain');
corrvals = arrayfun(@(x,y) corrmat(x,y), pairs(:,1), pairs(:,2));

switch distopt
    case 'distance'
        dist = ca_calc_interneuronal_distances(pairs, pos, 'distance');
    case 'channel'
        dist = ca_calc_interneuronal_distances(pairs, pos, 'channel');
end

%remove 0 dist
zeroidx = (dist == 0);
corrvals(zeroidx) = [];
dist(zeroidx) = [];

figure;
hold on
scatter(dist, corrvals);

X = [ones(length(dist),1) dist];
Y = corrvals;
[b, ~, ~,~, stats] = regress(Y,X);
x = min(dist)-5:max(dist)+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');

xlim([0, max(dist) + 50])

yl = ylim;
xl = xlim;

if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.8, yl(1) + (yl(2)-yl(1))*0.8, ...
        sprintf('R^{2} = %.2f\np < 0.001', stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.8, yl(1) + (yl(2)-yl(1))*0.8,...
        sprintf('R^{2} = %.2f\np = %.3f', stats(1), stats(3)))
end

xlabel('Pairwise distance (um)')
ylabel('Correlation values')

tickpref;    
temp = {stats corrvals dist};


for i = 1:nargout
    varargout{i} = temp{i};
end

