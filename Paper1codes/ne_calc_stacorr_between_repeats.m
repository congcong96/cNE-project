function corrvals = ne_calc_stacorr_between_repeats(exp_site_nedata, varargin)

stamat_main = ne_calc_sig_sta(exp_site_nedata);

for i = 1:length(varargin)
    
    stamat_rep = ne_calc_sig_sta(varargin{i}.exp_site_nedata);
    corrvals{i} = diag(corr(stamat_main',stamat_rep'));
end