function p = ne_batch_plot_NEmembers_depth(sig_sta)
% 
% NEmemdepths = cell(length(files),1);
% neudepths = cell(length(files),1);
% % NEdepths = cell(length(files),1);
% 
% 
% for i = 1:length(files)
%     load(files{i})
%     
%     [NEmemdepths{i}, neudepths{i}] = ne_calc_NEmembers_depth(exp_site_nedata);
% %     [NEdepths{i}] = ne_calc_NE_depth_median(exp_site_nedata);
%   
% end
% 
neurondepth = ne_calc_NE_statistics_by_group(sig_sta, 'depth');
    
% NEneurons = cell2mat(NEmemdepths);
% allneurons = cell2mat(neudepths);

% NEdepthsmat = cell2mat(NEdepths);
% NEneurons = cell2mat(NEmemdepths);
% allneurons = cell2mat(neudepths);

% maxval = max(allneurons);

% edges = 0:50:maxval+100;
% centers = (edges(1:end-1) + edges(2:end))/2;
% figure;
% hold on

% histogram(NEdepthsmat, edges, 'Normalization', 'probability');
% histogram(NEneurons, edges, 'Normalization', 'probability');
% histogram(allneurons, edges, 'Normalization', 'probability');

% NEdepthshist = histcounts(NEdepthsmat, edges, 'Normalization', 'probability');
% NEmemdepthshist = histcounts(NEneurons, edges, 'Normalization', 'probability');
% neudepthshist = histcounts(allneurons, edges, 'Normalization', 'probability');

% plot(centers, NEdepthshist);
% plot(centers, NEmemdepthshist, 'Color', color(1,:));
% plot(centers, neudepthshist, 'Color', color(2,:));
% legend(sprintf('depth of NE members, n = %d', length(NEneurons)),...
%     sprintf('depth of all recorded neurons, n = %d', length(allneurons)),...
%     'Location','Best')
% 
% xlabel('Depth (\mum)')
% ylabel('Probability')
% tickpref;
% print_mfilename(mfilename)

% color = eight_color_blind_palette('black', 'blue');
% 
% cat = [repmat({sprintf('All neurons, n = %d',length(allneurons))},...
%     length(allneurons), 1); repmat({sprintf('cNE member neurons, n = %d',...
%     length(NEneurons))}, length(NEneurons), 1)];
% g1 = gramm('x', [allneurons; NEneurons], 'color', cat);
% g1.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g1.set_names('x', 'Depth (\mum)', 'y', 'Ratio', 'color', 'Depth');
% g1.set_color_options('n_color', 2, 'n_lightness', 1, 'map', color);
% g1.set_order_options('x', 0);
% g1.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g1.draw();
% print_mfilename(mfilename);  
% 
% figure; hold on
% [f1, x1] = ecdf(allneurons);
% [f2, x2] = ecdf(NEneurons);
% plot(x1, f1, 'Color', color(1,:));
% plot(x2, f2, 'Color', color(2,:));
% legend(sprintf('depth of all recorded neurons, n = %d', length(allneurons)),...
%     sprintf('depth of NE members, n = %d', length(NEneurons)),...
%     'Location','Best')
% [~, p{1}] = kstest2(allneurons, NEneurons);



% color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');
% fn = fieldnames(neurondepth);
% 
% cat = cellfun(@(x, y) repmat({sprintf('%s, n = %d', x, length(y))},...
%     length(y), 1), fn, struct2cell(neurondepth), 'UniformOutput', 0);
% cat = vertcat(cat{:});
% g2 = gramm('x', cell2mat(struct2cell(neurondepth)), 'color', cat);
% g2.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g2.set_names('x', 'Depth (\mum)', 'y', 'Ratio', 'color', 'Depth');
% g2.set_color_options('n_color', 4, 'n_lightness', 1, 'map', color);
% g2.set_order_options('x', 0, 'color', 0);
% g2.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g2.draw();
% print_mfilename(mfilename);  

% 
% figure;
% hold on
% 
% for i = 1:length(fn)
%     depthhist = histcounts(neurondepth.(fn{i}), edges, 'Normalization', 'probability');
%     plot(centers, depthhist, 'Color', color(i,:));
% end
% legend(fn, 'Location', 'Best')
% 
% xlabel('Depth (\mum)')
% ylabel('Probability')
% tickpref;
% print_mfilename(mfilename)


% figure; hold on
% for i = 1:length(fn)
%     [f, x] = ecdf(neurondepth.(fn{i}));
%     plot(x, f, 'Color', color(i,:));
% end
% legend(fn, 'Location', 'Best');
% % legend(sprintf('All neurons, n = %d', length(allneurons)),...
% %     sprintf('cNE member neurons, n = %d', length(NEneurons)),...
% %     'Location', 'northeast')
% % [~, p] = kstest2(allneurons, NEneurons);
% % text(800, 0.4, sprintf('p = %.2f', p));
% tickpref;
% xlabel('Depth (\mum)')
% ylabel('Cumulative probability')
% 
% statmat = zeros(length(fn));
% comb = nchoosek(1:length(fn), 2);
% 
% for i = 1:size(comb,1)        
%     [~, tempp] = kstest2(neurondepth.(fn{comb(i,1)}), neurondepth.(fn{comb(i,2)}));
%     statmat(comb(i,1), comb(i,2)) = tempp * size(comb,1);
% end
% 
% statmat = statmat+statmat'; % to make statmat symmetrical    
% statmat(statmat > 1) = 1; % to make maximum value 1
% statmat(logical(eye(length(fn),length(fn)))) = nan; % to make diagonal nan
% 
% p{2} = statmat;    


color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');

figure;
axis ij
p = plot_plotspread_and_boxplot(neurondepth, 'kruskalwallis', color, -1);
ylabel('Depth (um)')
print_mfilename(mfilename);  

% axis ij


% histogram(NEdepthsmat,edges,'Normalization','probability');
% xlabel('Depth (\mum)')
% ylabel('Ratio')
% 
% x = xlim;
% y = ylim;
% 
% text(x(2)/4, y(2)/4*3, sprintf('n = %d', length(NEdepthsmat)))
% box on
% tickpref;
% hold off
% title('Depth of neurons in cNEs')
% print_mfilename(mfilename);
% 
% figure; hold on
% histogram(allneurons, edges, 'Normalization', 'probability');
% xlabel('Depth (\mum)')
% ylabel('Ratio')
% 
% text(x(2)/4, y(2)/4*3, sprintf('n = %d', length(allneurons)))
% box on
% tickpref;
% hold off
% title('Depth of all recorded neurons')
% print_mfilename(mfilename);
% 
% ylim([0 0.14])
% figure; hold on
% [fNE, xNE] = ecdf(NEdepthsmat);
% [fneu, xneu] = ecdf(allneurons);
% 
% plot(xNE, fNE, xneu, fneu)
% 
% legend('neurons in NEs','all recorded neurons','Location','best')
% xlabel('Depth (\mum)')
% ylabel('Cumulative probability')
% 
% [~,p] = kstest2(fNE,fneu);
% x = xlim;
% y = ylim;
% text(x(2)/4, y(2)/4*3, sprintf('p = %.2f', p))
% tickpref;
% title('CDF of depth of all neurons and cNE neurons')
% 
% hold off
% print_mfilename(mfilename);
