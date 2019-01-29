function STA_stats = ne_calc_STA_statistics_from_exp_site_nedata(exp_site_nedata)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;

neuronsta = nedata.stamat;
NEsta = nedata.NE_stamat;
nf = nedata.nf;
nlags = nedata.nlags;
spkcount = sum(nedata.spktrain, 2);
NEcount = sum(nedata.sta_NEtrain, 2);

% Calculate PTD
neuron_ptd = (max(neuronsta, [], 2) - min(neuronsta, [], 2)) ./ spkcount;
NE_ptd = (max(NEsta, [], 2) - min(NEsta, [], 2)) ./ NEcount;


% Calculate Moran's I
weightsmat = create_rf_spatial_autocorr_weights_matrix(nf, nlags, 1);

neuron_moransI = zeros(size(neuronsta, 1), 1);
for i = 1:length(neuron_moransI)
    neuron_moransI(i) = morans_i(reshape(neuronsta(i,:), nf, nlags), weightsmat);
end

NE_moransI = zeros(size(NEsta, 1), 1);
for i = 1:length(NE_moransI)
    NE_moransI(i) = morans_i(reshape(NEsta(i,:), nf, nlags), weightsmat);
end

% Calculate reliability index
stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);
neurontrain = nedata.sta_spktrain;
NEtrain = nedata.sta_NEtrain;
alltrain = [NEtrain; neurontrain];
% spkcount = sum(alltrain, 2);

[~, stabigmat, spkcountvec] = quick_calc_sta(stimstr.stimulus, alltrain, nedata.nlags, exp_site_nedata.stimlength);
rel_idx = calc_sta_reliability_index(stabigmat, spkcountvec); 

neuron_relidx = rel_idx(size(NEtrain,1)+1:end, :);
NE_relidx = rel_idx(1:size(NEtrain,1), :);

% Save to output struct array
STA_stats.neuron_ptd = neuron_ptd;
STA_stats.NE_ptd = NE_ptd;
STA_stats.neuron_moransI = neuron_moransI;
STA_stats.NE_moransI = NE_moransI;
STA_stats.neuron_reliability_idx = neuron_relidx;
STA_stats.NE_reliability_idx = NE_relidx;


end

