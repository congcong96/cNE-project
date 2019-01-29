function [ptd, moransI, rel_idx, sta] = ne_calc_STA_statistics_from_spktrain(spktrain, stim_mat, nlags, stimlength)

if ~exist('nlags','var')
    nlags = 20;
end

if ~exist('stimlength','var')
    stimlength = round(size(spktrain, 2) * 5 / 1000 / 60);
end

nf = size(stim_mat, 1);
sta = quick_calc_sta(stim_mat, spktrain, nlags);
spkcount = sum(spktrain, 2);

% Calculate PTD
ptd = (max(sta, [], 2) - min(sta, [], 2)) ./ spkcount;

% Calculate Moran's I
weightsmat = create_rf_spatial_autocorr_weights_matrix(nf, nlags, 1);

moransI = zeros(size(sta, 1), 1);
for i = 1:length(moransI)
    moransI(i) = morans_i(reshape(sta(i,:), nf, nlags), weightsmat);
end

% Calculate reliability index
[~, stabigmat, spkcountvec] = quick_calc_sta(stim_mat, spktrain, nlags, stimlength);
rel_idx = calc_sta_reliability_index(stabigmat, spkcountvec); 

end

