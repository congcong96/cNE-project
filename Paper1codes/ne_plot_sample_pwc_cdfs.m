function plot_sample_pwc_cdfs(pwcNE, pwcnNE, numbins)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%eliminate pwcs with meagre counts
count = {pwcNE.r12};
ridx = cellfun(@sum,count) < 200;
pwcNE(ridx) = [];
count = {pwcnNE.r12};
ridx = cellfun(@sum, count) < 200;
pwcnNE(ridx) = [];

if isempty(numbins)
    %get delays
    delay = pwcNE(1).delay;
    ddiff = delay(2)-delay(1);
    dbounds = [delay(1) delay(end)];

    %get number of bins to consider
    maxpd = max([abs([pwcNE.peakdelay]) abs([pwcnNE.peakdelay])]);
    numbins = (dbounds(2) - maxpd) / ddiff;
end

[NEad, NEcd, ~] = calc_pwc_sharpness_cdf(pwcNE, numbins, []);
[nNEad, nNEcd, ~] = calc_pwc_sharpness_cdf(pwcnNE, numbins, []);

for i = 1:5
    NEidx = randsample(size(NEcd,1), 1);
    nNEidx = randsample(size(nNEcd,1),1);
    
    figure;
    plot(NEad, NEcd(NEidx,:),'r', nNEad, nNEcd(nNEidx,:),'b');
    legend('neurons in same NE','neurons not in same NE')
    xlabel('Absolute delay - peak delay (ms)')
    ylabel('Cumulative probability')
    tickpref;
end
