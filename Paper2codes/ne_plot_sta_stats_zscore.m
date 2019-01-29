function ne_plot_sta_stats_zscore(sig_sta, plotopt)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

fn = fieldnames(sig_sta);
zscorefn = fn(contains(fn, 'zscore'));
NEzscorefn = zscorefn(contains(zscorefn, 'NE'));
neuronzscorefn = zscorefn(contains(zscorefn, 'neuron'));

sig_NE = [sig_sta.sig_NE];

for i = 1:length(NEzscorefn)
    temp = [sig_sta.(NEzscorefn{i})];    
    NEstruct.(sprintf('%s_sig', NEzscorefn{i})) = temp(sig_NE);
    NEstruct.(sprintf('%s_nonsig', NEzscorefn{i})) = temp(~sig_NE);
end

sig_neurons = cell2mat({sig_sta.sig_neurons}');

for i = 1:length(neuronzscorefn)
    temp = cell2mat({sig_sta.(neuronzscorefn{i})}');
    neuronstruct.(sprintf('%s_sig', neuronzscorefn{i})) = temp(sig_neurons);
    neuronstruct.(sprintf('%s_nonsig', neuronzscorefn{i})) = temp(~sig_neurons);
end



switch plotopt
    case 'xlog2d'
        %% cNE plot
        % get colormap
        cmin = min([NEstruct.NE_relidx_zscore_sig NEstruct.NE_relidx_zscore_nonsig]);
        cmax = max([NEstruct.NE_relidx_zscore_sig NEstruct.NE_relidx_zscore_nonsig]);
        cmap = disproportionate_colormap(cmin, cmax, 'cscheme', 'YlOrRd', 'trange', [0.2 1]);

        % plot transformed non-sig
        [tfx, tfy, tfz, reflectpt] = negsemilogx(NEstruct.NE_ptd_zscore_nonsig(:),...
            NEstruct.NE_moransI_zscore_nonsig(:), NEstruct.NE_relidx_zscore_nonsig(:),...
            'reflectpt', -1);
        figure('Position', [462, 252, 803, 581]);
        scatter(tfx, tfy, 12, tfz, 'o', 'filled');
        lennonsig = length(NEstruct.NE_ptd_zscore_nonsig);
        hold on

        % plot transformed sig
        scatter(log10(NEstruct.NE_ptd_zscore_sig) - reflectpt, NEstruct.NE_moransI_zscore_sig,...
            15,  NEstruct.NE_relidx_zscore_sig, '^', 'filled');
        lensig = length(NEstruct.NE_ptd_zscore_sig);

        % fix x axis
        xlab = abs(cellfun(@str2double, get(gca, 'XTickLabel')));
        xlabnew = xlab + reflectpt;
        yl = ylim;
        line([0 0], [yl(1) yl(2)], 'Color', 'k', 'LineStyle', '--')
        set(gca, 'XTickLabel', xlabnew);

        % standard plotting stuff
        colormap(cmap);
        xlabel('Log10(Peak-trough difference z-score)')
        ylabel('Moran''s I z-score')
        legend(sprintf('Non-significant cNEs, n = %d', lennonsig),...
            sprintf('Significant cNEs, n = %d', lensig),'Location', 'northwest')
        cb = colorbar;
        cb.Label.String = 'Reliability index z-score';
        cb.TickDirection = 'out';
        tickpref;

        %% Neuron plot
        % get colormap
        cmin = min([neuronstruct.neuron_relidx_zscore_sig; neuronstruct.neuron_relidx_zscore_nonsig]);
        cmax = max([neuronstruct.neuron_relidx_zscore_sig; neuronstruct.neuron_relidx_zscore_nonsig]);
        cmap = disproportionate_colormap(cmin, cmax, 'cscheme','YlOrRd', 'trange', [0.1 1]);

        % plot transformed non-sig
        [tfx, tfy, tfz, reflectpt] = negsemilogx(neuronstruct.neuron_ptd_zscore_nonsig(:),...
            neuronstruct.neuron_moransI_zscore_nonsig(:), neuronstruct.neuron_relidx_zscore_nonsig(:),...
            'reflectpt', -1);
        figure('Position', [462, 252, 803, 581]);
        scatter(tfx, tfy, 12, tfz, 'o', 'filled');
        lennonsig = length(neuronstruct.neuron_ptd_zscore_nonsig);
        hold on

        % plot transformed sig
        scatter(log10(neuronstruct.neuron_ptd_zscore_sig) - reflectpt, neuronstruct.neuron_moransI_zscore_sig,...
            15,  neuronstruct.neuron_relidx_zscore_sig, '^', 'filled');
        lensig = length(neuronstruct.neuron_ptd_zscore_sig);

        % fix x axis
        xlab = abs(cellfun(@str2double, get(gca, 'XTickLabel')));
        xlabnew = xlab + reflectpt;
        yl = ylim;
        line([0 0], [yl(1) yl(2)], 'Color', 'k', 'LineStyle', '--')
        set(gca, 'XTickLabel', xlabnew);

        % standard plotting stuff
        colormap(cmap);
        xlabel('Log10(Peak-trough difference z-score)')
        ylabel('Moran''s I z-score')
        legend(sprintf('Non-significant neurons, n = %d', lennonsig),...
             sprintf('Significant neurons, n = %d', lensig), 'Location', 'northwest')
        cb = colorbar;
        cb.Label.String = 'Reliability index z-score';
        cb.TickDirection = 'out';
        tickpref;

        
    case '2d'
        %% cNE plot
         % get colormap
        cmin = min([NEstruct.NE_relidx_zscore_sig NEstruct.NE_relidx_zscore_nonsig]);
        cmax = max([NEstruct.NE_relidx_zscore_sig NEstruct.NE_relidx_zscore_nonsig]);
        cmap = disproportionate_colormap(cmin, cmax, 'cscheme', 'YlOrRd', 'trange', [0.2 1]);

        % plot transformed non-sig
        figure('Position', [462, 252, 803, 581]);
        scatter(NEstruct.NE_ptd_zscore_nonsig, NEstruct.NE_moransI_zscore_nonsig,...
            12,  NEstruct.NE_relidx_zscore_nonsig, 'o', 'filled');
        lennonsig = length(NEstruct.NE_ptd_zscore_nonsig);
        hold on

        % plot transformed sig
        scatter(NEstruct.NE_ptd_zscore_sig, NEstruct.NE_moransI_zscore_sig,...
            15,  NEstruct.NE_relidx_zscore_sig, '^', 'filled');
        lensig = length(NEstruct.NE_ptd_zscore_sig);

        % standard plotting stuff
        colormap(cmap);
        xlabel('Log10(Peak-trough difference z-score)')
        ylabel('Moran''s I z-score')
        legend(sprintf('Non-significant cNEs, n = %d', lennonsig),...
            sprintf('Significant cNEs, n = %d', lensig),'Location', 'northwest')
        cb = colorbar;
        cb.Label.String = 'Reliability index z-score';
        cb.TickDirection = 'out';
        tickpref;
        
        %% Neuron plot
        cmin = min([neuronstruct.neuron_relidx_zscore_sig; neuronstruct.neuron_relidx_zscore_nonsig]);
        cmax = max([neuronstruct.neuron_relidx_zscore_sig; neuronstruct.neuron_relidx_zscore_nonsig]);
        cmap = disproportionate_colormap(cmin, cmax, 'cscheme','YlOrRd', 'trange', [0.1 1]);

        % plot transformed non-sig
        figure('Position', [462, 252, 803, 581]);
        scatter(neuronstruct.neuron_ptd_zscore_nonsig, neuronstruct.neuron_moransI_zscore_nonsig,...
            12,  neuronstruct.neuron_relidx_zscore_nonsig, 'o', 'filled');
        lennonsig = length(neuronstruct.neuron_ptd_zscore_nonsig);
        hold on

        % plot transformed sig
        scatter(neuronstruct.neuron_ptd_zscore_sig, neuronstruct.neuron_moransI_zscore_sig,...
            15,  neuronstruct.neuron_relidx_zscore_sig, '^', 'filled');
        lensig = length(neuronstruct.neuron_ptd_zscore_sig);

        % standard plotting stuff
        colormap(cmap);
        xlabel('Log10(Peak-trough difference z-score)')
        ylabel('Moran''s I z-score')
        legend(sprintf('Non-significant neurons, n = %d', lennonsig),...
             sprintf('Significant neurons, n = %d', lensig), 'Location', 'northwest')
        cb = colorbar;
        cb.Label.String = 'Reliability index z-score';
        cb.TickDirection = 'out';
        tickpref;       
                
        % to make it better spaced
        xlim([-20 140])
        
    case '3d'
        
        figure('Position', [462, 252, 803, 581]);
        scatter3(NEstruct.NE_ptd_zscore_nonsig, NEstruct.NE_moransI_zscore_nonsig,...
            NEstruct.NE_relidx_zscore_nonsig, 12, 'b', 'filled');
        lennonsig = length(NEstruct.NE_ptd_zscore_nonsig);
        
        hold on;        
        scatter3(NEstruct.NE_ptd_zscore_sig, NEstruct.NE_moransI_zscore_sig,...
            NEstruct.NE_relidx_zscore_sig, 12, 'r', 'filled');
        lensig = length(NEstruct.NE_ptd_zscore_sig);
        
        xlabel('Peak-trough difference z-score');
        ylabel('Moran''s I z-score');
        zlabel('Reliability index z-score');
        
        legend(sprintf('Non-significant cNE, n = %d', lennonsig),...
             sprintf('Significant cNE, n = %d', lensig), 'Location', 'southeast')
        
        figure('Position', [462, 252, 803, 581]);
        scatter3(neuronstruct.neuron_ptd_zscore_nonsig, neuronstruct.neuron_moransI_zscore_nonsig,...
            neuronstruct.neuron_relidx_zscore_nonsig, 12, 'b', 'filled');
        lennonsig = length(neuronstruct.neuron_ptd_zscore_nonsig);

        
        hold on;        
        scatter3(neuronstruct.neuron_ptd_zscore_sig, neuronstruct.neuron_moransI_zscore_sig,...
            neuronstruct.neuron_relidx_zscore_sig, 12, 'r', 'filled');
        lensig = length(neuronstruct.neuron_ptd_zscore_sig);

        
        xlabel('Peak-trough difference z-score');
        ylabel('Moran''s I z-score');
        zlabel('Reliability index z-score');
        xlim([-20 140])
        
        legend(sprintf('Non-significant neurons, n = %d', lennonsig),...
            sprintf('Significant neurons, n = %d', lensig), 'Location', 'southeast')
        
        
end
        
        
        
        