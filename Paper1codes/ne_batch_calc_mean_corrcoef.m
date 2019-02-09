function corrvec = ne_batch_calc_mean_corrcoef(files)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

corrvec = cell(length(files),1);

for i = 1:length(files)
    load(files{i}, 'exp_site_nedata')
    
    spktrain = exp_site_nedata.nedata.spktrain;
    corrmat = corr(spktrain');
    temp = triu(corrmat, 1);
    temp(temp == 0) = [];
    corrvec{i} = temp;
    
end
