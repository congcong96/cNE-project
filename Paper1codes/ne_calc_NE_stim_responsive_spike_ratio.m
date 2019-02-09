function respratio = ne_calc_NE_stim_responsive_spike_ratio(exp_site_nedata)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% if nargin == 1
%     plotopt = 0;
% end

nedata = exp_site_nedata.nedata;
NEmembers = nedata.NEmembers;
NEact = nedata.Activities;
spktrain = nedata.sta_spktrain;
NEthresh = nedata.NEthresh;

NEtrain = ne_upsample_NEact_using_member_neuron_activity(NEact, NEmembers, spktrain, NEthresh);

if ~isfield(exp_site_nedata, 'stimlength')
    stimlength = 10;
else
    stimlength = exp_site_nedata.stimlength;
end

if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
else
    dft = 10;
end

dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
nlags = nedata.nlags;
nf = nedata.nf;

[stimstr] = ca_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlength);

sta = quick_calc_sta(stimstr.stimulus, NEtrain, nlags, nf);

respratio = zeros(size(sta,1),1);

for i = 1:size(sta,1)

    [xprior, xposterior] = ne_sta_stimulus_projection(reshape(sta(i,:),nf, nlags), NEtrain(i,:), stimstr.stimulus);
    xprior(1:nlags-1) = [];
    [px,pspk,pxspk,xbinedges] = calc_px_pspk_pxspk(xprior,xposterior, 15);
    nl = pspk .* pxspk ./ px;
    xbins = edge2center(xbinedges);
    
%     [fiofit(i)] = hvv_fio_fit(xbins(:), nl, pspk);

    
    % find threshold of nonlinearity; this is the point where the nonlinearity
    % exceeds the average firing rate of the neuron (this uses the raw data, Brian's original)
    
    index = find(xbins > -1);
    bins_right = xbins(index);
    pspkx_right = nl(index);
    % bins_right = bins_right(1:end-2);
    % pspkx_right = pspkx_right(1:end-2);
    xi = linspace(min(bins_right), max(bins_right), 1000); % projection values
    yi = interp1(bins_right, pspkx_right, xi, 'spline'); % interpolated nonlinearity

    index = find(yi > pspk, 1);
    threshold = xi(index);
    
%     threshold = mean(threshold_temp);
        
    xmn = mean(xprior);
    xstd = std(xprior);
    xposterior_norm = (xposterior - xmn) ./ xstd;
    
    respratio(i) = sum(xposterior_norm > threshold) ./ length(xposterior_norm);
    
    
    
end
% 
% if plotopt == 1
%     
%     for i = 1:length(fiofit)
%         figure; hold on
%         plot(xbins, fiofit(i).pspkx, 'ko', 'markerfacecolor', 'k');
%         plot(fiofit(i).xFit, fiofit(i).pspkxFit, 'r-');
%         line([min(xbins) max(xbins)], [pspk pspk], 'Color', 'k', 'LineStyle', '--');
%         xmax = max([max(xbins) abs(min(xbins))]);
%         xmax = xmax + 2*xmax*0.05;
%         xlim([-xmax xmax]);
%         ylimit = get(gca,'ylim');
%         set(gca,'ylim', [-0.1*max(nl) max(ylimit)]);
%         set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%         legend('Data', 'HVV Fit', 'location', 'northwest');
%         title(sprintf('NMSE = %.3f, R2 = %.3f', fiofit(i).nmse, fiofit(i).r2));
%     end
% end   
    


