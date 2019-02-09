function ne_batch_plot_ICweights_fig3half(files)

figure;

minwt = zeros(length(files),1);
maxwt = zeros(length(files),1);

for i = 1:length(files)    
    load(files{i})
    
    pattern = exp_site_nedata.nedata.Patterns;
    NEmembers = exp_site_nedata.nedata.NEmembers;
    
    markers = cell(length(NEmembers),1);
    
    for j = 1:size(pattern,2)
        patsum = sum(pattern(:,j));
        markers{j} = [ones(length(NEmembers{j}),1)*j NEmembers{j}];

        if patsum < 0
            pattern(:,j) = -pattern(:,j);
        end
    end
    
    markersmat = cell2mat(markers);

    subplot(3,5,i)
    hold on
    imagesc(pattern);
    axis ij
    s = scatter(markersmat(:,1), markersmat(:,2), 10, [0 .8 0],'filled');
    s.Marker = 'square';
    
    box on
    set(gca, 'XTick', 1:size(pattern,2), 'XTickLabel', 1:size(pattern,2));
    tickpref;
    xlim([0.5 length(NEmembers)+0.5])
    ylim([0.5 size(pattern,1)+0.5])
    
    maxwt(i) = max(pattern(:));
    minwt(i) = min(pattern(:));
    
end

cbp = get(subplot(3,5,15),'Position');

cmapic = disproportionate_divergent_colormap(min(minwt),...
    max(maxwt), 1000, 'rdbu', 'brange', [0 0.7]);
colormap(cmapic);

colorbar('Position', [cbp(1)+cbp(3)+0.02  cbp(2)  0.02  cbp(2)+cbp(3)*5.7])

print_mfilename(mfilename);

