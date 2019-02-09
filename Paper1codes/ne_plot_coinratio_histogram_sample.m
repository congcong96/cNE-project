function ne_plot_coinratio_histogram_sample(exp_site_nedata, comb,numfiles, varargin)

fprintf('\n')

shufffiles = gfn('*shuffle*');
shufffiles = shufffiles(1:numfiles);

for i = 1:length(shufffiles)
    
    fprintf('Processing preISI file %d of %d...\n', i, length(shufffiles))
    
    load(shufffiles{i})
    spiketrains = surrspktrain.preISI;
    
    for j = 1:length(spiketrains)
        
        shuffcoinratio(i,j) = ne_calc_coincidence_within_spktrain(spiketrains{j}, comb);
        
    end
    
end

PRfiles = gfn('*-PR-*');
PRfiles = PRfiles(1:numfiles);

for i = 1:length(PRfiles)
    
    fprintf('Processing preRF file %d of %d...\n', i, length(PRfiles))
    
    load(PRfiles{i})
    spiketrains = surrspktrain.preRF;
    
    for j = 1:length(spiketrains)
        
        PRcoinratio(i,j) = ne_calc_coincidence_within_spktrain(spiketrains{j}, comb);
        
    end
    
end

DGfiles = gfn('*-ds10DG-*');
DGfiles = DGfiles(1:numfiles);

for k = 1:length(DGfiles)
    
    fprintf('Processing DG file %d of %d...\n', k, length(DGfiles))
    
    load(DGfiles{k})
    spiketrains = surrspktrain.DG;
    
    for ii = 1:length(spiketrains)
        
        DGcoinratio(k,ii) = ne_calc_coincidence_within_spktrain(spiketrains{ii}, comb);
        
    end
    
end

realcoinratio = ne_calc_coincidence_within_spktrain(exp_site_nedata.nedata.spktrain, comb);

[repmn] = ne_calc_coincidence_across_spktrains(exp_site_nedata.nedata.spktrain, varargin{1}, comb, 'fixedcomb'); %,1000,'median');

plotcolor = distinguishable_colors(5, 'w');

figure;
hold on

% edges = min(PPRcoinratio(:)): bindiff: max(PPRcoinratio(:)) + bindiff; 
h = histogram(shuffcoinratio, 'FaceColor', plotcolor(2,:), 'Normalization', 'probability');
shuffmn = median(shuffcoinratio(:));

bindiff = diff(h.BinEdges);

% edges = min(shuffcoinratio(:)): bindiff: max(shuffcoinratio(:)) + bindiff; 
% edges = min(repcoinratio(:)):bindiff:max(repcoinratio(:)) + bindiff;
% histogram(repcoinratio, edges, 'FaceColor', plotcolor(3,:), 'Normalization', 'probability')

edges = min(PRcoinratio(:)):bindiff:max(PRcoinratio(:)) + bindiff;
histogram(PRcoinratio, edges, 'FaceColor', plotcolor(4,:), 'Normalization', 'probability')
PRmn = median(PRcoinratio(:));


edges = min(DGcoinratio(:)):bindiff:max(DGcoinratio(:)) + bindiff;
histogram(DGcoinratio, edges, 'FaceColor', plotcolor(5,:), 'Normalization', 'probability')
DGmn = median(DGcoinratio(:));


y = ylim;

line([realcoinratio realcoinratio], [y(1) y(2)], 'Color', plotcolor(1,:), 'LineStyle', '--');

line([shuffmn shuffmn], [y(1) y(2)], 'Color', plotcolor(2,:), 'LineStyle', '--');

line([repmn repmn], [y(1) y(2)], 'Color', plotcolor(3,:), 'LineStyle', '--');

line([PRmn PRmn],[y(1) y(2)],'Color', plotcolor(4,:), 'LineStyle', '--');

line([DGmn DGmn], [y(1) y(2)], 'Color', plotcolor(5,:), 'LineStyle', '--');


xlabel('Coincidence ratio')
ylabel('Count')

legend('shuffle','PR', 'DG', 'real','shuffle (median)','repeat','PR (median)', 'DG (median)', 'Location','best')
box on
tickpref;

print_mfilename(mfilename);

hold off
fprintf('\n')
