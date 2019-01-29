function p = ne_compare_sig_vs_nonsig_cNE_spon_evoked_corr(evokedfiles, sponfiles)

assert(length(evokedfiles) == length(sponfiles))

sig_corrvals = cell(length(evokedfiles), 1);
nonsig_corrvals = cell(length(evokedfiles), 1);

for i = 1:length(evokedfiles)
    
    load(evokedfiles{i})
    
    evokedpat = exp_site_nedata.nedata.Patterns;
    NEsig = exp_site_nedata.nedata.sig_NE_sta;
    
    load(sponfiles{i})
    
    sponpat = exp_site_nedata.nedata.Patterns;
    
    NE_sim = abs(corr(evokedpat, sponpat));
    corrvals = max(NE_sim, [], 2);
    
    sig_corrvals{i} = corrvals(NEsig);
    nonsig_corrvals{i} = corrvals(~NEsig);
    

end

datastruct.sig = cell2mat(sig_corrvals);
datastruct.nonsig = cell2mat(nonsig_corrvals);
figure;
p = plot_plotspread_and_boxplot(datastruct);

