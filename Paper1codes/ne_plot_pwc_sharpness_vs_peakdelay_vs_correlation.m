function varargout = ne_plot_pwc_sharpness_vs_peakdelay_vs_correlation(exp_site_nedata, spktrain, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


fsd = 2000;
stimdur = 1/fsd*1000 * size(spktrain,2);

if nargin == 3
    pwc = varargin{1};
elseif nargin < 3
    pwc = pairwisecorr_function(spktrain, stimdur, fsd);
else
    error('Too many inputs')
end

[~, ~, interpvals, pwcidx] = calc_pwc_sharpness_cdf(pwc, 40, 0.5, 0);
interpvals(isnan(interpvals)) = 0;

newpwc = pwc(pwcidx);

peakdelay = abs([newpwc.peakdelay]');
pairs = cell2mat({newpwc.pairs}');
corrmat = corr(exp_site_nedata.nedata.spktrain');

corrvals = arrayfun(@(x,y) corrmat(x,y),pairs(:,1),pairs(:,2));


% pwcpeak = zeros(length(newpwc),1);
% for k = 1:length(newpwc)    
%     pwcpeak(k) = max(newpwc(k).r12 ./ sum(newpwc(k).r12));    
% end


graphtitle = sprintf('%s-site%d-%s',exp_site_nedata.exp, exp_site_nedata.site, exp_site_nedata.stim);
%correct underscores
graphtitle = regexprep(graphtitle, '_', '\\_');

figure;
set(gcf, 'Position', [200 50 1500 900])

subplot(3,1,1)
posidx = corrvals >= 0;
% [f, gof{1}] = fit(log(corrvals(posidx)),interpvals(posidx),'exp1');
% plot(f,log(corrvals(posidx)),interpvals(posidx))

[f, gof{1}] = fit(corrvals,interpvals,'exp1');
plot(f,corrvals,interpvals)

xlabel('spike train corr coeff')
ylabel('CC sharpness')
title(graphtitle)

xlim([-.05 .4])
x = xlim;
y = ylim;
text((x(2) - x(1))/4*3+x(1), (y(2) - y(1)) /5*3+y(1), sprintf('R^{2} = %.3f', gof{1}.adjrsquare));  
tickpref;

lgd = legend;
lgd.Location = 'best';
box off


subplot(3,1,2)
% [f, gof{2}] = fit(log(corrvals(posidx)),peakdelay(posidx),'exp1');
% plot(f,log(corrvals(posidx)),peakdelay(posidx))
[f, gof{2}] = fit(corrvals,peakdelay,'exp1');
plot(f,corrvals,peakdelay)

xlabel('spike train corr coeff')
ylabel('CC peak delay')
title(graphtitle)

xlim([-.05 .4])
x = xlim;
y = ylim;
text((x(2) - x(1))/4*3+x(1), (y(2) - y(1)) /5*3+y(1), sprintf('R^{2} = %.3f', gof{2}.adjrsquare));  
tickpref;

lgd = legend;
lgd.Location = 'best';
box off

subplot(3,1,3)
[f, gof{3}] = fit(interpvals, peakdelay, 'exp1');
plot(f, interpvals, peakdelay)

xlabel('CC sharpness')
ylabel('CC peak delay')
title(graphtitle)

x = xlim;
y = ylim;
text((x(2) - x(1))/6+x(1), (y(2) - y(1)) /5*3+y(1), sprintf('R^{2} = %.3f', gof{3}.adjrsquare));  
tickpref;

lgd = legend;
lgd.Location = 'best';
box off


% subplot(3,2,4)
% [f, gof{4}] = fit(log(pwcpeak), interpvals, 'poly1');
% plot(f, log(pwcpeak), interpvals)
% 
% xlabel('log (PWC peak value)')
% ylabel('PWC sharpness')
% title(graphtitle)
% 
% x = xlim;
% y = ylim;
% text((x(2) - x(1))/5+x(1), (y(2) - y(1)) /4+y(1), sprintf('R^{2} = %.3f', gof{4}.adjrsquare));  
% tickpref;
% 
% lgd = legend;
% lgd.Location = 'best';
% 
% 
% subplot(3,2,5)
% [f, gof{5}] = fit((pwcpeak), peakdelay, 'exp1');
% plot(f, (pwcpeak), peakdelay)
% 
% xlabel('log (PWC peak value)')
% ylabel('PWC peak delay')
% title(graphtitle)

% x = xlim;
% y = ylim;
% text((x(2) - x(1))/5+x(1), (y(2) - y(1)) /4*3+y(1), sprintf('R^{2} = %.3f', gof{5}.adjrsquare));  
% tickpref;
% 
% lgd = legend;
% lgd.Location = 'best';

% subplot(3,2,6)
% [f, gof{6}] = fit(log(pwcpeak(posidx)), log(corrvals(posidx)), 'poly1');
% plot(f, log(pwcpeak(posidx)), log(corrvals(posidx)))
% 
% xlabel('log (PWC peak value)')
% ylabel('spiketrain dot product')
% title(graphtitle)
% 
% x = xlim;
% y = ylim;
% text((x(2) - x(1))/5+x(1), (y(2) - y(1)) /4*3+y(1), sprintf('R^{2} = %.3f', gof{6}.adjrsquare));  
% tickpref;
% 
% lgd = legend;
% lgd.Location = 'best';
% X = [ones(length(corrvals),1) corrvals];
% Y = interpval;
% [b, ~, ~,~, stats] = regress(Y,X);
% x = min(corrvals)-.05:0.001:max(corrvals)+.05;
% y = b(1) + b(2)*x;
% plot(x,y,'k--');





% yl = ylim;
% xl = xlim;

% if stats(3) < 0.001
%     text(xl(1)+(xl(2)-xl(1))*0.2, yl(1) + (yl(2)-yl(1))*0.2, ...
%         sprintf('R^{2} = %.2f\np < 0.001', stats(1)))
% else
%     text(xl(1)+(xl(2)-xl(1))*0.2, yl(1) + (yl(2)-yl(1))*0.2,...
%         sprintf('R^{2} = %.2f\np = %.3f', stats(1), stats(3)))
% end
% tickpref;    
temp = {gof corrvals interpvals peakdelay pwc};

print_mfilename(mfilename);

for i = 1:nargout
    varargout{i} = temp{i};
end

