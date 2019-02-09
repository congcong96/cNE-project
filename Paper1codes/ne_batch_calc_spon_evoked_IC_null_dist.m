function corrvalstruct = ne_batch_calc_spon_evoked_IC_null_dist(sponfiles, evokedfiles)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

assert(length(sponfiles) == length(evokedfiles))

corrvalstruct(length(sponfiles)).sponfile = [];

for i = 1:length(sponfiles)
    load(sponfiles{i}, 'exp_site_nedata')
    spon_nedata = exp_site_nedata.nedata;
    load(evokedfiles{i}, 'exp_site_nedata')
    evoked_nedata = exp_site_nedata.nedata;
    clear exp_site_nedata
    
    [shuffledcorrvals, sponcorrvals, evokedcorrvals] = ne_calc_spon_evoked_IC_null_dist_thresholded(spon_nedata, evoked_nedata);
    
    corrvalstruct(i).sponfile = sponfiles{i};
    corrvalstruct(i).evokedfile = evokedfiles{i};
    corrvalstruct(i).shuffledcorrvals = shuffledcorrvals;
    corrvalstruct(i).sponcorrvals = sponcorrvals;
    corrvalstruct(i).evokedcorrvals = evokedcorrvals;

end

