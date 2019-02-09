function [strfhw, cchw] = ne_plot_strfcorr_vs_pwc_examples(strf, trigger, pwc, plotopt)
% strf_overlap STRF overlap for two neurons
% 
%     pairs_plot_strf_overlap_crosscorr_example(strf, trigger, ccpairs)
% 
%     strf : struct array of strfs, where some channels have two STRFs.
%
%     trigger : vector of trigger times, in sample number.
%
%     ccpairs : struct array holding pairwise correlation data for the
%     channels in STRF. ccpairs holds the cross-covariance function, the
%     correlogram, the peak delay, half-width, etc. from correlation
%     analysis.
%
%     
%     The function finds channels where there are at least two STRFs, and
%     then correlates the STRFs using a slice through the CF.
%
%     The RF slice correlation is then compared to the spike train correlation
if nargin == 3
    plotopt = 0;
end

faxis = strf(1).faxis;
taxis = strf(1).taxis;
dtemp = taxis(2) - taxis(1);
if ( dtemp < 0.1 )
    taxis = taxis * 1000; % convert to ms
end

dt = ( taxis(2) - taxis(1) ); % size of bins
maxlag = ceil(20 / dt); % Number of bins to cover 50 ms



% Earlier code:
% % Part 1. Get slice through peaks:
% 
% [i1, j1] =  find(max(max(strf1)) == strf1);
% slice1 = strf1(i1,:);
% slice11 = slice1 ./ max([abs(slice1) eps]); % strf1, peak1
% 
% [i2, j2] =  find(max(max(strf2)) == strf2);
% slice2 = strf2(i2,:);
% slice22 = slice2 ./ max([abs(slice2) eps]); % strf2, peak 2
% 
% 
% % slice2 = strf2(i1,:);
% % slice21 = slice2 ./ max([abs(slice2) eps]); % strf 2, peak 1
% % 
% % slice1 = strf1(i2,:)
% % slice12 = slice1 ./ max([abs(slice1) eps]); % strf 1, peak 2
% 
% % [slice11(:) slice22(:) slice21(:) slice12(:)]
% 
% 
% % [c12, lags] = xcorr(slice11, slice22, maxlag, 'coeff');
% % [c11, lags] = xcorr(slice11, slice21, maxlag, 'coeff');
% % [c22, lags] = xcorr(slice12, slice22, maxlag, 'coeff');
% 
% [c12, lags] = xcorr(slice11, slice22, maxlag);
% % [c11, lags] = xcorr(slice11, slice21, maxlag);
% % [c22, lags] = xcorr(slice12, slice22, maxlag);
% c12 = c12 ./ max(c12);
% % c11 = c11 ./ max(c11);
% % c22 = c22 ./ max(c22);
% 
% 
% lags = lags * dt;
% 
% 
% index12 = find(c12 == max(c12));
% peak12 = c12(index12);
% lag_peak12 = lags(index12);
% 
% index_low12 = find(c12 < 0.5 *peak12 & lags < lag_peak12 );
% index_low12 = max([1 index_low12]);
% lag_low12 = lags( max( index_low12 ) );
% 
% index_high12 = find(c12 < 0.5 *peak12 & lags > lag_peak12 );
% index_high12 = min([length(c12) index_high12]);
% lag_high12 = lags( min( index_high12 ) );
% 
% 
% 
% 
% clf;
% 
% % Part 1.
% subplot(4,3,1);
% plot_strf_symmetric_colormap(strf1);
% title('STRF1');
% 
% subplot(4,3,2);
% plot_strf_symmetric_colormap(strf2);
% title('STRF2');
% 
% subplot(4,3,4);
% plot(taxis, slice11, 'k-');
% xlim([min(taxis) max(taxis)]);
% ylim([min(slice11) max(slice11)]);
% box on;
% tickpref;
% title('STRF1 peak, STRF1 slice');
% 
% subplot(4,3,5);
% plot(taxis, slice22, 'k-');
% xlim([min(taxis) max(taxis)]);
% ylim([min(slice22) max(slice22)]);
% box on;
% tickpref;
% title('STRF2 peak, STRF2 slice');
% 
% subplot(4,3,6);
% hold on;
% patch([lag_low12 lag_low12 lag_high12 lag_high12],1*[-1 1 1 -1],0.85*[1 1 1]);
% patch(1.25*[-1 -1 1 1],1*[-1 1 1 -1],[0.6 0.6 0.6]);
% plot(lags,c12,'k-');
% xtick = sort([-50 -25 -10 10 25 50]);
% set(gca,'xtick', xtick, 'xticklabel', xtick);
% xlim([min(lags) max(lags)]);
% box on;
% tickpref;
% title('Slice through each peak');





p = 0.01;


% chan = [strf.chan];
% chan_unique = unique(chan);


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


% for i = 1:length(chan_unique)
% 
% chan_unique(i)
% 
%   index = find( chan_unique(i) == chan );
% 
%   cmb = nchoosek(index, 2); % determine all possible pairwise combinations
% 
%   [nr, nc] = size(cmb);

%   for j = 1:nr

strfhw = zeros(size(pairs,1),1);
cchw = zeros(size(pairs,1),1);

for j = 1:size(pairs,1)
    
    index1 = pairs(j,1);
    index2 = pairs(j,2);

%     depth = strf(index1).depth;
%     position = strf(index1).position;
%     stim = strf(index1).stim;
%     atten = strf(index1).atten;
% 
%     fs = strf(index1).fs;


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
    taxispos = taxis(indt);
    rfsig1 = rfsig1(:,indt);
    rfsig2 = rfsig2(:,indt);

    indf = find(faxis <= 10000);
    faxispos = faxis(indf);
    rfsig1 = rfsig1(indf,:);
    rfsig2 = rfsig2(indf,:);


    % Get the cross-cov function for the pair of neurons
%     cc = pairs_get_crosscorr_for_strf_pair(strf(index1), strf(index2), ccpairs);
    cc = pwclist(j);
    cc_delay = cc.delay;
    cc_q12 = cc.q12;
    cc_q12 = cc_q12 ./ max(cc_q12);
    

%     % Correlate the receptive fields
%     [i1, j1] =  find(max(max(rfsig1)) == rfsig1);
% 
% 
%     slice1 = rfsig1(i1,:);

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

%     index12 = find(c12 == max(c12));
%     peak12 = c12(index12);
%     lag_peak12 = lags(index12);
% 
%     index_low12 = find(c12 < 0.5 *peak12 & lags < lag_peak12 );
%     index_low12 = max([1 index_low12]);
%     lag_low12 = lags( max( index_low12 ) );
% 
%     index_high12 = find(c12 < 0.5 *peak12 & lags > lag_peak12 );
%     index_high12 = min([length(c12) index_high12]);
%     lag_high12 = lags( min( index_high12 ) );
%     
    strfhw(j,:) = find_strf_corr_halfwidth(lags,c12);
    
    ccfit = fit(cc_delay(:), cc_q12(:), 'gauss1');
%     cchw(j,:) = find_spike_xcorr_fit_halfwidth(cc_delay, ccfit);

    hb = bar(cc_delay,cc_q12);
    set(hb, 'facecolor', 'k');
    xtick = sort([-50 -25 -10 0 10 25 50]);
    set(gca,'xtick', xtick, 'xticklabel', xtick);
    xlim([min(lags) max(lags)]);
    ylim(1.03*[ min(min([c12(:); cc_q12(:)])) max(max([c12(:); cc_q12(:)]))]);
    box on;
    cchw(j,:) = find_crude_spike_xcorr_halfwidth(cc_q12(lagstartidx:lagendidx), dt, 2);
    
    if plotopt == 1

        clf;

        % Part 1.
        subplot(3,2,1);
        hold on;
        plot_strf_symmetric_colormap(rfsig1); %,taxispos, faxispos/1000);
%         plot(xlim, [i1 i1], 'k-');
        xlim([0 length(taxispos)]);
        ylim([0 length(faxispos)]);
        set(gca,...
        'xtick', [1 round(length(taxispos)/2) length(taxispos)], ...
        'xticklabel', [0 50 100], ...
        'ytick', [1 round(length(faxispos)/2) length(faxispos)], ...
        'yticklabel', [5 10 20]);
        title('STRF1');
        ylabel('Frequency (kHz)');
        box on;
        title(sprintf('Neuron %d',index1))

        subplot(3,2,3);
        hold on;
        plot_strf_symmetric_colormap(rfsig2); %,taxispos, faxispos/1000);
%         plot(xlim, [i2 i2], 'k-');
        xlim([0 length(taxispos)]);
        ylim([0 length(faxispos)]);
        set(gca,...
        'xtick', [1 round(length(taxispos)/2) length(taxispos)], ...
        'xticklabel', [0 50 100], ...
        'ytick', [1 round(length(faxispos)/2) length(faxispos)], ...
        'yticklabel', [5 10 20]);
        title('STRF2');
        xlabel('Time (ms)');
        ylabel('Frequency (kHz)');
        title(sprintf('Neuron %d',index2))
        box on;

        subplot(3,2,2);
        plot(taxispos, slice11, 'k-');
        xlim([min(taxispos) max(taxispos)]);
        ylim([min(slice11) max(slice11)]);
        set(gca,...
        'xtick', [min(taxispos) max(taxispos)/2 max(taxispos)], ...
        'xticklabel', [0 50 100]);
        box on;
        tickpref;
        title('STRF1 mean');

        subplot(3,2,4);
        plot(taxispos, slice22, 'k-');
        xlim([min(taxispos) max(taxispos)]);
        ylim([min(slice22) max(slice22)]);
        set(gca,...
        'xtick', [min(taxispos) max(taxispos)/2 max(taxispos)], ...
        'xticklabel', [0 50 100]);
        box on;
        tickpref;
        xlabel('Time (ms)');
        title('STRF2 mean');

%         subplot(3,2,5);
%         hold on;
%         patch([lag_low12 lag_low12 lag_high12 lag_high12],1*[-1 1 1 -1],0.85*[1 1 1]);
%         patch(1.25*[-1 -1 1 1],1*[-1 1 1 -1],[0.6 0.6 0.6]);
%         plot(lags,c12,'k-');
%         xtick = sort([-50 -25 -10 0 10 25 50]);
%         set(gca,'xtick', xtick, 'xticklabel', xtick);
%         xlim([min(lags) max(lags)]);
%         box on;
%         tickpref;
%         title('Slice through each peak');

        subplot(3,2,5);
        plot(ccfit, cc_delay, cc_q12)
%         ccfit = loess(cc_delay,cc_q12,min(cc_delay):0.1:max(cc_delay),0.25, 1);
%         ccfit = ccfit./max(ccfit);
%         hold on
%         scatter(cc_delay, cc_q12, 6, 'r', '.');
%         plot(min(cc_delay):0.1:max(cc_delay), ccfit);
        ylim([0 1])
        xlim([min(lags) max(lags)]);        


        subplot(3,2,6);
        hold on;
        plot(lags,c12,'k-');
        hb = bar(cc_delay,cc_q12);
        set(hb, 'facecolor', 'k');
        xtick = sort([-50 -25 -10 0 10 25 50]);
        set(gca,'xtick', xtick, 'xticklabel', xtick);
        xlim([min(lags) max(lags)]);
        ylim(1.03*[ min(min([c12(:); cc_q12(:)])) max(max([c12(:); cc_q12(:)]))]);
        box on;
        tickpref;
    %     legend('STRF Overlap', 'Spike Train Corr');
        xlabel('Delay (ms)');

        set(gcf,'position', [417 200 1000 600]);

        pause
    end

% % end % (for j)

end % (for i)


return;








