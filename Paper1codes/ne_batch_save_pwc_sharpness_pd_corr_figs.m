function batch_ne_save_pwc_sharpness_pd_corr_figs(NEfiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

base = regexp(NEfiles, '^\S+(?=(-ne-))','match','once');
cmbfiles = cellfun(@(x) [x '.mat'],base,'UniformOutput',0);

for i = 1:length(NEfiles)
    load(NEfiles{i})
    load(cmbfiles{i})
    
    if exist('allpwc','var')
        ne_plot_pwc_sharpness_vs_correlation(exp_site_nedata, spktrain,allpwc);
    else
        [~, ~, ~, ~, allpwc] = ne_plot_pwc_sharpness_vs_correlation(exp_site_nedata, spktrain);
    end
    
    try
        save(cmbfiles{i}, 'edges','trigger','spk','position','spktrain','strf','allpwc')
    catch
        save(cmbfiles{i}, 'edges','trigger','spk','position','spktrain','allpwc')
    end
    
    savefig(gcf, sprintf('I:\\Cell_Assemblies\\PaperNEs\\corrval_pwcsharpness_pwcpd\\%s',base{i}));
    
    clearvars -except NEfiles base cmbfiles
    close all

end

