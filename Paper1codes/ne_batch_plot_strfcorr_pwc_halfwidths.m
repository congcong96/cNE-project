function ne_batch_plot_strfcorr_pwc_halfwidths(nefiles, spkfiles)

assert(length(nefiles) == length(spkfiles));
for i = 1:length(nefiles)
    
    fprintf('\nProcessing file %d of %d',i,length(nefiles))
    nebasefile = regexp(nefiles{i}, '^\S+(?=(-fs))','match','once');
    spkbasefile = regexp(spkfiles{i}, '^\S+(?=(-fs))','match','once');
    assert(strcmp(nebasefile, spkbasefile));
    
    load(nefiles{i})
    load(spkfiles{i})
    
    hw = ne_calc_strfcorr_and_pwc_halfwidths(strf, trigger, allpwc, exp_site_nedata);
    ne_plot_strfcorr_vs_pwc_halfwidth_comparison(hw);

    
    
end
