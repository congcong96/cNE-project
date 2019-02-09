function [NEhw, nNEhw, p] = ne_batch_calc_NE_vs_nonNE_strf_halfwidth(strffiles)

p = zeros(length(strffiles),1);

figure;

for i = 1:length(strffiles)
    
    fprintf('Processing %s...\n', strffiles{i})
    
    load(strffiles{i}, 'strf', 'trigger')
    basefilename = regexp(strffiles{i}, '^\S+(?=(.mat))', 'match', 'once');
    nefile = [basefilename '-ne-20dft.mat'];
    
    load(nefile, 'exp_site_nedata')
    subplot(4,4,i)
    [NEhw{i}, nNEhw{i}, p(i)] = ne_calc_NE_vs_nonNE_strf_halfwidth(exp_site_nedata, strf, trigger, 1);
    
end
print_mfilename(mfilename);
fprintf('\n')