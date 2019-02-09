function p = ne_batch_plot_NE_depth_median_and_span(sig_sta)

% NEdepths = cell(length(files),1);
% NEspans = cell(length(files),1);
% 
% alldepths = cell(length(files),1);
% allspans = cell(length(files),1);
% 
% for i = 1:length(files)
%     load(files{i})
%     
%     [NEdepths{i}, NEspans{i}] = ne_calc_NE_depth_median(exp_site_nedata);
%     
%     sitedepths = cell2mat(exp_site_nedata.nedata.position');
%     sitedepths = sitedepths(:,2);
%     alldepths{i} = repmat(median(sitedepths), length(NEdepths{i}), 1);
%     allspans{i} = repmat(max(sitedepths) - min(sitedepths), length(NEspans{i}), 1);
%   
% end

% NEdepthsmat = cell2mat(NEdepths);
% NEspansmat = cell2mat(NEspans);

% alldepthsmat = cell2mat(alldepths);
% allspansmat = cell2mat(allspans);

% edges = min(alldepthsmat):50:max(alldepthsmat)+50;
% wedges = min(allspansmat):50:max(allspansmat)+50;

% color = eight_color_blind_palette('black', 'blue');

% %% PDF median depth all neurons vs cNE
% cat = [repmat({'All'}, length(alldepthsmat), 1); repmat({'cNEs'}, ...
%     length(NEdepthsmat), 1)];
% g1 = gramm('x', [alldepthsmat; NEdepthsmat], 'color', cat);
% g1.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g1.set_names('x', 'Depth (\mum)', 'y', 'Ratio', 'color', 'Depth');
% g1.set_color_options('n_color', 2, 'n_lightness', 1, 'map', color);
% g1.set_order_options('x', 0);
% g1.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g1.draw();
% print_mfilename(mfilename);  
% 
% %% CDF median depth all neurons vs cNE
% figure; hold on
% [f1, x1] = ecdf(alldepthsmat);
% [f2, x2] = ecdf(NEdepthsmat);
% plot(x1, f1, 'Color', color(1,:));
% plot(x2, f2, 'Color', color(2,:));
% legend({'Median depth of recorded neurons','Median depth of cNE members'},...
%     'Location','Best')
% [~, p{1}] = kstest2(allspansmat, NEspansmat);
% 
% %% PDF span all neurons vs cNE
% cat = [repmat({'All'}, length(allspansmat), 1); repmat({'cNEs'}, ...
%     length(NEspansmat), 1)];
% g2 = gramm('x', [allspansmat; NEspansmat], 'color', cat);
% g2.stat_bin('edges', wedges, 'normalization', 'probability', 'geom', 'line');
% g2.set_names('x', 'Span (\mum)', 'y', 'Ratio', 'color', 'Span');
% g2.set_color_options('n_color', 2, 'n_lightness', 1, 'map', color);
% g2.set_order_options('x', 0);
% g2.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g2.draw();
% print_mfilename(mfilename);  
% 
% %% CDF span all neurons vs cNE
% figure; hold on
% [f1, x1] = ecdf(allspansmat);
% [f2, x2] = ecdf(NEspansmat);
% plot(x1, f1, 'Color', color(1,:));
% plot(x2, f2, 'Color', color(2,:));
% legend({'Span of recorded neurons','Span of cNEs'},...
%     'Location','Best')
% [~, p{2}] = kstest2(allspansmat, NEspansmat);
% 
%% PDF median depth all groups
medianstruct = ne_calc_NE_statistics_by_group(sig_sta, 'median');

color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');
% fn = fieldnames(medianstruct);
% 
% cat = cellfun(@(x, y) repmat({sprintf('%s, n = %d', x, length(y))},...
%     length(y), 1), fn, struct2cell(medianstruct), 'UniformOutput', 0);
% cat = vertcat(cat{:});
% g3 = gramm('x', cell2mat(struct2cell(medianstruct)), 'color', cat);
% g3.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g3.set_names('x', 'Depth (\mum)', 'y', 'Ratio', 'color', 'Depth');
% g3.set_color_options('n_color', 4, 'n_lightness', 1, 'map', color);
% g3.set_order_options('x', 0, 'color', 0);
% g3.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g3.draw();
% print_mfilename(mfilename);
% 
% %% CDF median depth all groups
% figure; hold on
% for i = 1:length(fn)
%     [f, x] = ecdf(medianstruct.(fn{i}));
%     plot(x, f, 'Color', color(i,:));
% end
% legend(fn, 'Location', 'Best');
% tickpref;
% xlabel('Depth (\mum)')
% ylabel('Cumulative probability')
% 
% statmat = zeros(length(fn));
% comb = nchoosek(1:length(fn), 2);
% 
% for i = 1:size(comb,1)        
%     [~, tempp] = kstest2(medianstruct.(fn{comb(i,1)}), medianstruct.(fn{comb(i,2)}));
%     statmat(comb(i,1), comb(i,2)) = tempp * size(comb,1);
% end
% 
% statmat = statmat+statmat'; % to make statmat symmetrical    
% statmat(statmat > 1) = 1; % to make maximum value 1
% statmat(logical(eye(length(fn),length(fn)))) = nan; % to make diagonal nan
% 
% p{3} = statmat; 

%% Box plot median depth all groups
figure;
axis ij
p{1} = plot_plotspread_and_boxplot(medianstruct, 'kruskalwallis', color, -1);
print_mfilename(mfilename);
ylabel('cNE depth (um)')


%% PDF span all groups
spanstruct = ne_calc_NE_statistics_by_group(sig_sta, 'span');
% 
% color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');
% fn = fieldnames(spanstruct);
% 
% cat = cellfun(@(x, y) repmat({sprintf('%s, n = %d', x, length(y))},...
%     length(y), 1), fn, struct2cell(spanstruct), 'UniformOutput', 0);
% cat = vertcat(cat{:});
% g4 = gramm('x', cell2mat(struct2cell(spanstruct)), 'color', cat);
% g4.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g4.set_names('x', 'Span (\mum)', 'y', 'Ratio', 'color', 'Span');
% g4.set_color_options('n_color', 4, 'n_lightness', 1, 'map', color);
% g4.set_order_options('x', 0, 'color', 0);
% g4.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g4.draw();
% print_mfilename(mfilename);
% 
% %% CDF span all groups
% figure; hold on
% for i = 1:length(fn)
%     [f, x] = ecdf(spanstruct.(fn{i}));
%     plot(x, f, 'Color', color(i,:));
% end
% legend(fn, 'Location', 'Best');
% tickpref;
% xlabel('Span (\mum)')
% ylabel('Cumulative probability')
% 
% statmat = zeros(length(fn));
% comb = nchoosek(1:length(fn), 2);
% 
% for i = 1:size(comb,1)        
%     [~, tempp] = kstest2(spanstruct.(fn{comb(i,1)}), spanstruct.(fn{comb(i,2)}));
%     statmat(comb(i,1), comb(i,2)) = tempp * size(comb,1);
% end
% 
% statmat = statmat+statmat'; % to make statmat symmetrical    
% statmat(statmat > 1) = 1; % to make maximum value 1
% statmat(logical(eye(length(fn),length(fn)))) = nan; % to make diagonal nan
% 
% p{5} = statmat; 

%% Box plot span all groups
figure;
p{2} = plot_plotspread_and_boxplot(spanstruct, 'kruskalwallis', color, 1);
print_mfilename(mfilename);
ylabel('cNE span (um)')

% figure;
% hold on
% 
% histogram(NEdepthsmat,edges,'Normalization','probability');
% % histogram(categorical(alldepthsmat),'Normalization','probability');
% xlabel('Depth (\mum)')
% ylabel('Ratio')
% x = xlim;
% y = ylim;
% text(x(2)/4, y(2)/4*3, sprintf('n = %d', length(NEdepthsmat)))
% box on
% tickpref;
% hold off
% print_mfilename(mfilename);
% 
% figure;
% hold on
% histogram(NEspansmat,wedges,'Normalization','probability');
% xlabel('Span (\mum)')
% ylabel('Ratio')
% 
% x = xlim;
% y = ylim;
% 
% text(x(2)/4*3, y(2)/4*3, sprintf('n = %d', length(NEspansmat)))
% box on
% tickpref;
% print_mfilename(mfilename);
