function p = ne_batch_plot_distance_histogram(sig_sta)

% NEdist = cell(length(files),1);
% alldist = cell(length(files),1);
% 
% for i = 1:length(files)
%     load(files{i})
% %     basename = regexp(files{i}, '^\S+(?=(-ne-))','match','once');
% %     spkfile = [basename '.mat'];
% %     load(spkfile, 'spk')
%     
%     [NEdist{i}, alldist{i}] = ne_plot_distance_histogram(exp_site_nedata, 'distance', 0);   
%   
% end

% NEdistmat = cell2mat(NEdist);
% alldistmat = cell2mat(alldist);

diststruct = ne_calc_NE_statistics_by_group(sig_sta, 'pairwise_distance');



% switch distopt
%     case 'distance'
% edges = 0:50:max(NEdistmat);
%         histogram(NEdistmat,edges);
%         xlabel('Pairwise distance (\mum)')
%         ylabel('Count')
%         title('Distance histogram of neurons within NEs')
%         tickpref;
%         hold off
%         
%         figure;
%         histogram(NEdistmat,edges,'Normalization','probability');
%         xlabel('Pairwise distance (\mum)')
%         ylabel('Ratio')
%         x = xlim;
%         y = ylim;
%         text(x(2)/2, y(2)/2, sprintf('n = %d', length(NEdistmat)))
%         tickpref;
%         print_mfilename(mfilename);
%         
%         figure;
%         histogram(alldistmat,edges,'Normalization','probability');
%         xlabel('Pairwise distance (\mum)')
%         ylabel('Ratio')
%         x = xlim;
%         ylim([y(1) y(2)])
%         text(x(2)/2, y(2)/2, sprintf('n = %d', length(alldistmat)))
%         tickpref;
%         print_mfilename(mfilename);
% color = eight_color_blind_palette('black', 'blue');
% 
% cat = [repmat({sprintf('All pairs, n = %d',length(alldistmat))},...
%     length(alldistmat), 1); repmat({sprintf('cNE pairs, n = %d',...
%     length(NEdistmat))}, length(NEdistmat), 1)];
% g = gramm('x', [alldistmat; NEdistmat], 'color', cat);
% g.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g.set_names('x', 'Pairwise distance', 'y', 'Ratio', 'color', 'Pairs');
% g.set_color_options('n_color', 2, 'n_lightness', 1, 'map', color);
% g.set_order_options('x', 0, 'color', 0);
% g.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g.draw();
% print_mfilename(mfilename);     
% 
% figure; hold on
% [f1, x1] = ecdf(alldistmat);
% [f2, x2] = ecdf(NEdistmat);
% plot(x1, f1, 'Color', color(1,:));
% plot(x2, f2, 'Color', color(2,:));
% legend(sprintf('All pairs, n = %d', length(alldistmat)),...
%     sprintf('cNE pairs, n = %d', length(NEdistmat)),...
%     'Location', 'northeast')
% [~, p{1}] = kstest2(alldistmat, NEdistmat);
% text(800, 0.4, sprintf('p = %.3e', p{1}));
% tickpref;
% xlabel('Pairwise distance (\mum)')
% ylabel('Cumulative probability')
% 
% 
% color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');
% fn = fieldnames(diststruct);
% 
% cat = cellfun(@(x, y) repmat({sprintf('%s, n = %d', x, length(y))},...
%     length(y), 1), fn, struct2cell(diststruct), 'UniformOutput', 0);
% cat = vertcat(cat{:});
% g2 = gramm('x', cell2mat(struct2cell(diststruct)), 'color', cat);
% g2.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g2.set_names('x', 'Pairwise distance', 'y', 'Ratio', 'color', 'Depth');
% g2.set_color_options('n_color', 4, 'n_lightness', 1, 'map', color);
% g2.set_order_options('x', 0, 'color', 0);
% g2.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g2.draw();
% print_mfilename(mfilename);  
% 
% 
% figure; hold on
% for i = 1:length(fn)
%     [f, x] = ecdf(diststruct.(fn{i}));
%     plot(x, f, 'Color', color(i,:));
% end
% legend(fn, 'Location', 'Best');
% % legend(sprintf('All neurons, n = %d', length(allneurons)),...
% %     sprintf('cNE member neurons, n = %d', length(NEneurons)),...
% %     'Location', 'northeast')
% % [~, p] = kstest2(allneurons, NEneurons);
% % text(800, 0.4, sprintf('p = %.2f', p));
% tickpref;
% xlabel('Pairwise distance (\mum)')
% ylabel('Cumulative probability')
% 
% statmat = zeros(length(fn));
% comb = nchoosek(1:length(fn), 2);
% 
% for i = 1:size(comb,1)        
%     [~, tempp] = kstest2(diststruct.(fn{comb(i,1)}), diststruct.(fn{comb(i,2)}));
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
p = plot_plotspread_and_boxplot(diststruct, 'kruskalwallis', color, 1);
ylabel('Pairwise distance (um)')
print_mfilename(mfilename);  


%     case 'channel'
%         NEdistmat = categorical(NEdistmat);
%         histogram(NEdistmat)
%         xlabel('Number of channels apart')
%         ylabel('Count')
%         title('Distance histogram of neurons within NEs')
%         tickpref;
%         hold off
%         
%         figure;
%         h = histogram(NEdistmat,'Normalization','probability');
%         xlabel('Number of channels apart')
%         ylabel('Ratio')
% 
%         x = xlim;
%         y = ylim;
% 
%         text(x(2)/2, y(2)/2, sprintf('n = %d', length(NEdistmat)))
%         tickpref;
%         print_mfilename(mfilename);
% end


