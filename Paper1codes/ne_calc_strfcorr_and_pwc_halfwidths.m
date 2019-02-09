function hw = ne_calc_strfcorr_and_pwc_halfwidths(strf, trigger, pwc, exp_site_nedata)

faxis = strf(1).faxis;
taxis = strf(1).taxis;
dtemp = taxis(2) - taxis(1);
if ( dtemp < 0.1 )
    taxis = taxis * 1000; % convert to ms
end

dt = ( taxis(2) - taxis(1) ); % size of bins
maxlag = ceil(20 / dt); % Number of bins to cover 50 ms

p = 0.01;

pwclist = [];

if iscell(pwc)
    for i = 1:length(pwc)
        pwclist = [pwclist pwc{i}];
    end
    pairs = cell2mat({pwclist.pairs}');
%     crosscorrhw = cell2mat({pwclist.halfwidth}');
    
else
    pairs = cell2mat({pwc.pairs}');
%     crosscorrhw = [pwc.halfwidth];
    pwclist = pwc;
    
end

counts = cellfun(@sum, {pwclist.r12});
removeidx = counts < 200;
pwclist(removeidx) = [];
pairs(removeidx, :) = [];


% find NE and non-NE pairs
NEneurons = ne_find_NE_pairs_or_groups(exp_site_nedata, 2);

[~, NEidx] = intersect(pairs,NEneurons,'rows');
nonNEidx = setdiff(1:length(pwclist), NEidx)';

hw.NE(length(NEidx)).strfhw = [];
hw.nonNE(length(nonNEidx)).strfhw = [];

nidx = 1;
nnidx = 1;

for j = 1:size(pairs,1)
    
    index1 = pairs(j,1);
    index2 = pairs(j,2);

    % Now get the STRF similarity data
    rf1 = strf(index1).rfcontra;
    rf2 = strf(index2).rfcontra;

    n01 = strf(index1).n0contra;
    n02 = strf(index2).n0contra;

    mdb = strf(index1).mdb;

    fs = strf(index1).fs;

    stimdur = ( trigger(end) - trigger(1) ) / fs;

    soundtype = strf(index1).stim;

    rfsig1 = significant_strf(rf1, p, n01, mdb, stimdur, soundtype);
    rfsig2 = significant_strf(rf2, p, n02, mdb, stimdur, soundtype);

%     rfsig1 = rf1;
%     rfsig2 = rf2;

    indt = find(taxis >= -5 & taxis <= 95);
    rfsig1 = rfsig1(:,indt);
    rfsig2 = rfsig2(:,indt);

    indf = find(faxis <= 10000);
    rfsig1 = rfsig1(indf,:);
    rfsig2 = rfsig2(indf,:);


    % Get the cross-cov function for the pair of neurons
%     cc = pairs_get_crosscorr_for_strf_pair(strf(index1), strf(index2), ccpairs);
    cc = pwclist(j);
    cc_delay = cc.delay;
    cc_q12 = cc.q12;
    cc_q12 = cc_q12 ./ max(cc_q12);


    % find mean of RF and normalize
    slice1 = mean(rfsig1);
    slice11 = slice1 ./ max([abs(slice1) eps]); % strf1, peak1

%     [i2, j2] =  find(max(max(rfsig2)) == rfsig2);
%     slice2 = rfsig2(i2,:);

    slice2 = mean(rfsig2);
    slice22 = slice2 ./ max([abs(slice2) eps]); % strf2, peak 2

    [c12, lags] = xcorr(slice11, slice22, maxlag);
    c12 = c12 ./ max(c12);

    lags = lags * dt;
    
    lagstartidx = find(min(lags) == cc_delay);
    lagendidx = find(max(lags) == cc_delay);
    
%     hb = bar(cc_delay,cc_q12);
%     set(hb, 'facecolor', 'k');
%     xtick = sort([-50 -25 -10 0 10 25 50]);
%     set(gca,'xtick', xtick, 'xticklabel', xtick);
%     xlim([min(lags) max(lags)]);
%     ylim(1.03*[ min(min([c12(:); cc_q12(:)])) max(max([c12(:); cc_q12(:)]))]);
%     box on;


    if any(j == NEidx)
        
        hw.NE(nidx).strfhw = find_strf_corr_halfwidth(lags,c12);
        hw.NE(nidx).cchw = find_crude_spike_xcorr_halfwidth(cc_q12(lagstartidx:lagendidx), dt, 6);
        nidx = nidx + 1;
        
    elseif any(j == nonNEidx)
        
        hw.nonNE(nnidx).strfhw = find_strf_corr_halfwidth(lags,c12);
        hw.nonNE(nnidx).cchw = find_crude_spike_xcorr_halfwidth(cc_q12(lagstartidx:lagendidx), dt, 6);
        nnidx = nnidx + 1;
        
    else
        
        keyboard
        
    end
    
    
end

return