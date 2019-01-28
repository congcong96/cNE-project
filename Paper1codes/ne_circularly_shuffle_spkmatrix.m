function surrspkmat = ne_circularly_shuffle_spkmatrix (exp_site_nedata, opt)

% Circularly shuffle neuronal spiketrain

% Updated 7/20/16 by js

if ~exist('opt','var')
    opt = 'regular';
end

switch opt
    case 'regular'
        spktrain = exp_site_nedata.nedata.spktrain;
    case 'sta'
        spktrain = exp_site_nedata.nedata.sta_spktrain;
end

surrspkmat = zeros(size(spktrain));

for i = 1:size(spktrain,1)
    surrspkmat(i,:) = circularly_shuffle_neuronal_spktrain (spktrain(i,:));
end

