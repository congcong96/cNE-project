function p = ne_plot_spon_evoked_IC_corr_zscore(corrvalstruct, sigopt)

sig = cell(length(corrvalstruct), 1);
nonsig = cell(length(corrvalstruct), 1);

for i = 1:length(corrvalstruct)
    
    shuffcorrvals = corrvalstruct(i).shuffledcorrvals;
    
    if sigopt
        thresh = 0.01;
        evokedcorrvals = corrvalstruct(i).evokedcorrvals;
        pval = zeros(length(evokedcorrvals), 1);
        for j = 1:length(evokedcorrvals)
            compvec = [abs(shuffcorrvals); abs(evokedcorrvals(j))];
            pval(j) = sum(abs(evokedcorrvals(j)) <= compvec) ./ length(compvec);
        end
        idx = pval < thresh;
    else
        idx = true(length(corrvalstruct(i).evokedcorrvals), 1);
    end
    
    miu = mean(shuffcorrvals);
    sigma = std(shuffcorrvals);
    normevoked = (corrvalstruct(i).evokedcorrvals - miu) ./ sigma;
    
    load(corrvalstruct(i).evokedfile, 'exp_site_nedata')
    sig_NE = exp_site_nedata.nedata.sig_NE_sta;
    
    sig{i} = abs(normevoked(sig_NE & idx));
    nonsig{i} = abs(normevoked(~sig_NE & idx));
    
end

datastruct.sig = cell2mat(sig')';
datastruct.nonsig = cell2mat(nonsig')';

color = eight_color_blind_palette('orange','blue');
figure;
p = plot_plotspread_and_boxplot(datastruct, 'ranksum', color, 1);
ylabel('Spontaneous-evoked z-scored correlation score')
    