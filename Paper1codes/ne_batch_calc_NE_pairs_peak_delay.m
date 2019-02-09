function pd = ne_batch_calc_NE_pairs_peak_delay(nefiles)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pd = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    load(nefiles{i}, 'exp_site_nedata')
    pwcNE = exp_site_nedata.nedata.pwc;
    pwcNE = cell2mat(pwcNE);
    
    pwcNE([pwcNE.significant] == 0) = [];
    pd{i} = [pwcNE.peakdelay];

end

