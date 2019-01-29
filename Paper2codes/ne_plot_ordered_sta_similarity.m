function ne_plot_ordered_sta_similarity(exp_site_nedata)

nedata = exp_site_nedata.nedata;
stamat = nedata.stamat;
NEmembers = nedata.NEmembers;

NEmemvec = cell2mat(NEmembers);
NEmem_stamat = stamat(NEmemvec, :);
simmat = corr(NEmem_stamat');
simmat = set_matrix_diagonal(simmat, 0);

%zero other 1s
simmat(simmat > 0.9999 & simmat < 1.0001) = 0;

NEcumsize = [0; cumsum(cellfun('length', NEmembers))] + 0.5;
cmap = brewermap(1000, 'Reds');


figure; hold on
imagesc(simmat);
colormap(cmap)

for i = 1:length(NEcumsize) - 1
    
    line([NEcumsize(i) NEcumsize(i+1) NEcumsize(i+1) NEcumsize(i) NEcumsize(i)],...
        [NEcumsize(i) NEcumsize(i) NEcumsize(i+1) NEcumsize(i+1) NEcumsize(i)],...
        'Color', 'k', 'LineWidth', 2)

end

axis ij
hold off

xlim([0.5 size(simmat,1) + 0.5])
ylim([0.5 size(simmat,1) + 0.5])