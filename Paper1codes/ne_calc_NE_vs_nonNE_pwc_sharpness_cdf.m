function [hdNE, hdnNE] = ne_calc_NE_vs_nonNE_pwc_sharpness_cdf(exp_site_nedata, allpwc, pval, plotopt, numbins)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    pval = 0.5; 
    plotopt = 0;
    numbins = [];
elseif nargin == 3
    plotopt = 0;
    numbins = [];
elseif nargin == 4
    numbins = [];
end

pwcNE = exp_site_nedata.nedata.pwc;
pwcNE = cell2mat(pwcNE);

NEpairs = cell2mat({pwcNE.pairs}');
[~, uniqueidx] = unique(NEpairs, 'rows', 'stable');

pwcNE = pwcNE(uniqueidx);

% nonNEpairs = ne_find_non_NE_pairs_or_groups(exp_site_nedata, 2, length(pwcNE)*1.5);
allpwcpairs = cell2mat({allpwc.pairs}');
[~, nonNEidx] = setdiff(allpwcpairs, NEpairs, 'rows');

pwcnNE = allpwc(nonNEidx);

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

%get delay value at pval
[~, ~, hdNE] = calc_pwc_sharpness_cdf(pwcNE, numbins, pval);
[~, ~, hdnNE] = calc_pwc_sharpness_cdf(pwcnNE, numbins, pval);

hdNE(isnan(hdNE)) = 0;
hdnNE(isnan(hdnNE)) = 0;

if plotopt == 1
    
    ne_plot_NE_vs_nonNE_pwc_cdf_interp(hdNE, hdnNE);
    
end

end

