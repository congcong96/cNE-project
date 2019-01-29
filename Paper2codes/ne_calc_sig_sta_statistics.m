function [sta_sig, ptd, moransI, rel_idx, siglevel] = ne_calc_sig_sta_statistics(sta, locator, stim, nreps, nlags, stimlength, pval)

% rtfopt, varargin)
% rand_test_sta_locator Significance test for STA calculated using a locator
% 
%    [sta_sig, siglevel, rand_dist] = ...
%       rand_test_sta_locator(sta, locator, stimulus, nreps, pval)
%
%    computes a randomization test on STA to determine the significant pixels
%    in STA. The test is performed by shifting the spike train, in LOCATOR, by
%    an amount equal to 
% 
%       shiftsize = round( length(locator)/(nreps+1) );
% 
%    and then calculating a random STA. NREPS determines the shift size of the
%    spike train, and also the size of random distribution. PVAL is the
%    significance level of the test.
% 
%    STA_SIG is the significant part of STA at the level PVAL. SIGLEVEL is 
%    the level that was used to determine signifcant parts of STA, and
%    RAND_DIST is the distribution of randomized STA values. It is a vector 
%    of length = NREPS*length(STA(:))
% 
%    Inputs
%    ----------------------------------------------------------
%    sta : spike-triggered average previously estimated
%    locator : spike train of 1's and 0's
%    stim : matrix of the stimulus. Each row is one trial.
%    nreps : number randomization to perform
%    pval : significance level
% 
%    Outputs
%    ----------------------------------------------------------
%    sta_sig : significant parts of sta
%    siglevel : value used to threshold sta
%    rand_dist : distribution used to obtain siglevel
% 
%    Adapted from ne_sig_sta_from_stim_obs_resp.m by on 6/18/18 by JS.

if ~exist('pval', 'var')
    pval = 95;
end

nf = size(stim, 1);

% get fixed shift size
shiftsize = round( length(locator)/(nreps+1) );

stalength = size(sta,2); % nf x nlags
% dim: num_units x (nf x nlags) x nreps
sta_rand_mat = zeros(size(sta,1),stalength, nreps); 
% peak-trough difference (one value per unit per shuffled iteration)
ptd = zeros(size(locator,1), nreps);
% reliability index(100 values per unit per shuffled iteration)
rel_idx = cell(1, nreps);
% morans I (one value per unit per shuffled iteration)
moransI = zeros(size(locator,1), nreps);

spkcount = sum(locator, 2);

% initialize moran's I weights matrix
weightsmat = create_rf_spatial_autocorr_weights_matrix(nf, nlags, 1);
     
fprintf('\nCalculating significant STAs...\n')

for i = 1:nreps
    
    fprintf('\nIteration %d of %d', i, nreps)
    % circularly shift spike train
    shift = i * shiftsize;
    loc_rand = circshift(locator, [0 shift]);
    
    % calculate one-minute chunks of STA (and overall STA)
    [temp, stabigmat, spkcountvec] = quick_calc_sta(stim, loc_rand, nlags, stimlength);
    % calculate reliability index
    rel_idx{i} = calc_sta_reliability_index(stabigmat, spkcountvec); 
    % calculate raw peak-trough difference
    ptd(:,i) = (max(temp, [], 2) - min(temp, [], 2)) ./ spkcount;
    
    % calculate Moran's I
    for j = 1:size(locator,1)
        moransI(j, i) = morans_i(reshape(temp(j,:), nf, nlags), weightsmat);
    end

    sta_rand_mat(:,:,i) = temp;
    
end

rel_idx = cell2mat(rel_idx);

lower = (100 - pval)/2;
upper = 100 - lower;

% get significance threshold for each STA
siglevel = prctile(reshape(sta_rand_mat, [size(sta_rand_mat,1), ...
    size(sta_rand_mat, 2) * size(sta_rand_mat, 3)]),[lower upper],2);

sta_sig = zeros(size(sta));
% contigcount = cell(size(sta,1), 1);

for j = 1:size(sta,1)
    % get thresholded STA
    sta_sig(j,  sta(j,:) < siglevel (j,1) | sta(j,:) > siglevel(j,2)) = ...
        sta(j, sta(j,:) < siglevel (j,1) | sta(j,:) > siglevel(j,2));
    
end

fprintf('\n')
return;






