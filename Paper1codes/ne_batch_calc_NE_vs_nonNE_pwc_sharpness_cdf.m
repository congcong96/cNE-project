function [allhdNE, allhdnNE] = ne_batch_calc_NE_vs_nonNE_pwc_sharpness_cdf(files)

allhdNE = cell(length(files),1);
allhdnNE = cell(length(files),1);

for i = 1:length(files)
    fprintf('\nProcessing %s...\n', files{i})
    load(files{i}, 'exp_site_nedata')
    filebase = regexp(files{i}, '^\S+(?=(-ne))', 'match', 'once');
    
    spkfile = [filebase '.mat'];
    load(spkfile, 'allpwc')
    
    [allhdNE{i}, allhdnNE{i}] = ne_calc_NE_vs_nonNE_pwc_sharpness_cdf(exp_site_nedata, allpwc, 0.5, 0, 20);
    
    clear('exp_site_nedata','allpwc')
    
end