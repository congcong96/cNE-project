function allhw = ne_batch_calc_strfcorr_and_pwc_halfwidths(nefiles, spkfiles)

assert(length(nefiles) == length(spkfiles));

allhw.NE = [];
allhw.nonNE = [];
for i = 1:length(nefiles)
    
    fprintf('\nProcessing file %d of %d',i,length(nefiles))
    nebasefile = regexp(nefiles{i}, '^\S+(?=(-fs))','match','once');
    spkbasefile = regexp(spkfiles{i}, '^\S+(?=(-fs))','match','once');
    assert(strcmp(nebasefile, spkbasefile));
    
    load(nefiles{i})
    load(spkfiles{i})
    
    hw = ne_calc_strfcorr_and_pwc_halfwidths(strf, trigger, allpwc, exp_site_nedata);
    allhw.NE = [allhw.NE hw.NE];
    allhw.nonNE = [allhw.nonNE hw.nonNE];
    
end

fprintf('\n')
   
    