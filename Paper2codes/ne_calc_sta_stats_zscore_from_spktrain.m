function stastats = ne_calc_sta_stats_zscore_from_spktrain(spktrain, stim_mat, nreps, fs)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('nreps','var')
    nreps = 100;
end

if ~exist('fs','var')
    fs = 200; % 5 ms: default STA calculation resolution
end

stimlen = round(size(spktrain,2)/fs/60);

nlags = 20;
nf = size(stim_mat, 1);

%% Calculate real spike train metrics

[sta, sta_chunks, spkcountvec] = quick_calc_sta(stim_mat, spktrain, nlags, stimlen);
spkcount = sum(spktrain, 2);

stastats.spkcount = spkcount;
% Calculate PTD
stastats.ptd = (max(sta, [], 2) - min(sta, [], 2)) ./ spkcount;

% Calculate Moran's I
weightsmat = create_rf_spatial_autocorr_weights_matrix(nf, nlags, 1);

moransI = zeros(size(sta, 1), 1);
for i = 1:length(moransI)
    moransI(i) = morans_i(reshape(sta(i,:), nf, nlags), weightsmat);
end
stastats.moransI = moransI;

% Calculate reliability_index
stastats.reliability_idx = calc_sta_reliability_index(sta_chunks, spkcountvec); 


%% Calculate null distribution of metrics

% get fixed shift size
shiftsize = round(size(spktrain, 2)/(nreps+1));

% peak-trough difference (one value per unit per shuffled iteration)
shuff_ptd = zeros(size(spktrain,1), nreps);

% reliability index(100 values per unit per shuffled iteration)
shuff_reliability_idx = cell(1, nreps);

% morans I (one value per unit per shuffled iteration)
shuff_moransI = zeros(size(spktrain,1), nreps);

fprintf('\nCalculating null distributions...\n')

for i = 1:nreps
    
    fprintf('\nIteration %d of %d', i, nreps)
    % circularly shift spike train
    shift = i * shiftsize;
    loc_rand = circshift(spktrain, [0 shift]);
    
    % calculate one-minute chunks of STA (and overall STA)
    [sta, sta_chunks, spkcountvec] = quick_calc_sta(stim_mat, loc_rand, nlags, stimlen);
    % calculate reliability index
    shuff_reliability_idx{i} = calc_sta_reliability_index(sta_chunks, spkcountvec); 
    % calculate raw peak-trough difference
    shuff_ptd(:,i) = (max(sta, [], 2) - min(sta, [], 2)) ./ spkcount;
    
    % calculate Moran's I
    for j = 1:size(spktrain,1)
        shuff_moransI(j, i) = morans_i(reshape(sta(j,:), nf, nlags), weightsmat);
    end
    
end

stastats.shuff_ptd = shuff_ptd;
stastats.shuff_moransI = shuff_moransI;
shuff_reliability_idx = cell2mat(shuff_reliability_idx);
stastats.shuff_reliability_idx = shuff_reliability_idx;


%% Calculate zscores of real data

ptd_mu = mean(shuff_ptd, 2, 'omitnan');
ptd_sigma = std(shuff_ptd, 0, 2, 'omitnan');
stastats.ptd_zscore = (stastats.ptd - ptd_mu)./ ptd_sigma;

reliability_idx_mu = mean(shuff_reliability_idx, 2, 'omitnan');
reliability_idx_sigma = std(shuff_reliability_idx, 0, 2, 'omitnan');
stastats.reliability_idx_zscore = (mean(stastats.reliability_idx, 2, 'omitnan') - reliability_idx_mu)./ reliability_idx_sigma;

moransI_mu = mean(shuff_moransI, 2, 'omitnan');
moransI_sigma = std(shuff_moransI, 0, 2, 'omitnan');
stastats.moransI_zscore = (stastats.moransI - moransI_mu) ./ moransI_sigma;

%% Calculate pvals for real data

ptd_pval = zeros(size(shuff_ptd, 1), 1);
reliability_idx_pval = zeros(size(shuff_ptd, 1), 1);
moransI_pval = zeros(size(shuff_ptd, 1), 1);

for i = 1:size(shuff_ptd, 1)
    
    ptd_pval(i) = (sum(stastats.ptd(i) <= shuff_ptd(i,:)) + 1) ./ (length(shuff_ptd(i,:)) + 1);
    reliability_idx_pval(i) = (sum(stastats.reliability_idx(i) <= shuff_reliability_idx(i,:)) + 1)...
        ./ (length(shuff_reliability_idx(i,:)) + 1);
    moransI_pval(i) = (sum(stastats.moransI(i) <= shuff_moransI(i,:)) + 1) ./ (length(shuff_moransI(i,:)) + 1);

end

stastats.ptd_pval = ptd_pval;
stastats.reliability_idx_pval = reliability_idx_pval;
stastats.moransI_pval = moransI_pval;

