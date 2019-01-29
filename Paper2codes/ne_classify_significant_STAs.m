function [neuron_sig, NE_sig, details] = ne_classify_significant_STAs(exp_site_nedata, varargin)

ip = inputParser;
addRequired(ip, 'exp_site_nedata', @isstruct);
addParameter(ip, 'plotopt', [], @(x) isempty(x) || strcmp(x, 'raw') || strcmp(x, 'sig'));
addParameter(ip, 'pvalthresh', [], @isstruct);
% addParameter(ip, 'ptdthresh', 0.05, @(x) x <= 1.1 || x >= 0);
% addParameter(ip, 'moransIthresh', 0.05, @(x) x <= 1.1 || x >= 0);
% addParameter(ip, 'relidxthresh', 0.05, @(x) x <= 1.1 || x >= 0);
addParameter(ip, 'zscorecheck', 0, @(x) x == 0 || x == 1);
parse(ip, exp_site_nedata, varargin{:});

exp_site_nedata = ip.Results.exp_site_nedata;
plotopt = ip.Results.plotopt;
% ptdthresh = ip.Results.ptdthresh;
% moransIthresh = ip.Results.moransIthresh;
% relidxthresh = ip.Results.relidxthresh;
pvalthresh = ip.Results.pvalthresh;
zscorecheck = ip.Results.zscorecheck;

if isempty(pvalthresh)
    pvalthresh.ptd_thresh = 0.05;
    pvalthresh.moransI_thresh = 0.05;
    pvalthresh.reliability_idx_thresh = 0.05;
else
    assert(pvalthresh.ptd_thresh <= 1.1 && pvalthresh.ptd_thresh >= 0)
    assert(pvalthresh.moransI_thresh <= 1.1 && pvalthresh.moransI_thresh >= 0)
    assert(pvalthresh.reliability_idx_thresh <= 1.1 && pvalthresh.reliability_idx_thresh >= 0)
end


nedata = exp_site_nedata.nedata;

zscorethresh = 4;

fn = fieldnames(nedata.shuffled_stats);

for i = 1:length(fn)
    
    nulldist = nedata.shuffled_stats.(fn{i});
    realval = mean(nedata.STA_stats.(fn{i}), 2, 'omitnan');
    mu = mean(nulldist, 2, 'omitnan');
    sigma = std(nulldist, 0, 2, 'omitnan');
    z_score = (realval - mu) ./ sigma;
    
    dist = [nulldist realval];
    stattype = regexp(fn{i}, '(?<=^\w{2,6}_)\S+$', 'match', 'once');
    pval = sum(realval <= dist, 2) ./ size(dist, 2);
    sig.(fn{i}) = pval < pvalthresh.(sprintf('%s_thresh', stattype));
    sigz.(fn{i}) = z_score >= zscorethresh;
    
    details.([fn{i} '_zscore']) = z_score;
    details.([fn{i} '_pval']) = pval;
    
end
    
% 
% % compare PTD of real filters and shuffled filters
% NE_ptddist = nedata.shuffled_stats.NE_ptd;
% NE_ptd = nedata.STA_stats.NE_ptd;
% NE_ptd_mu = mean(NE_ptddist, 2, 'omitnan');
% NE_ptd_sigma = mean(NE_ptddist, 2
% 
% NE_ptddist = [NE_ptddist NE_ptd];
% NE_ptd_pval = sum(NE_ptd <= NE_ptddist, 2) ./ size(NE_ptddist, 2);
% NE_ptd_sig = NE_ptd_pval < ptdthresh;
% NE_ptd_zscore = nedata.sig_classifying_stats.NE_ptd_zscore;
% NE_ptd_sigz = NE_ptd_zscore >= zscorethresh;
% 
% neuron_ptddist = nedata.shuffled_stats.neuron_ptd;
% neuron_ptd = nedata.STA_stats.neuron_ptd;
% neuron_ptddist = [neuron_ptddist neuron_ptd];
% neuron_ptd_pval = sum(neuron_ptd <= neuron_ptddist, 2) ./ size(neuron_ptddist, 2);
% neuron_ptd_sig = neuron_ptd_pval < ptdthresh;
% neuron_ptd_zscore = nedata.sig_classifying_stats.neuron_ptd_zscore;
% neuron_ptd_sigz = neuron_ptd_zscore >= zscorethresh;
% 
% % compare Moran's I real and shuffled filters
% NE_midist = nedata.shuffled_stats.NE_moransI;
% NE_mi = nedata.STA_stats.NE_moransI;
% NE_midist = [NE_midist NE_mi];
% NE_mi_pval = sum(NE_mi <= NE_midist, 2) ./ size(NE_midist, 2);
% NE_mi_sig = NE_mi_pval < moransIthresh;
% NE_mi_zscore = nedata.sig_classifying_stats.NE_moransI_zscore;
% NE_mi_sigz = NE_mi_zscore >= zscorethresh;
% 
% neuron_midist = nedata.shuffled_stats.neuron_moransI;
% neuron_mi = nedata.STA_stats.neuron_moransI;
% neuron_midist = [neuron_midist neuron_mi];
% neuron_mi_pval = sum(neuron_mi <= neuron_midist, 2) ./ size(neuron_midist, 2);
% neuron_mi_sig = neuron_mi_pval < moransIthresh;
% neuron_mi_zscore = nedata.sig_classifying_stats.neuron_moransI_zscore;
% neuron_mi_sigz = neuron_mi_zscore >= zscorethresh;
% 
% % compare reliability index of real and shuffled filters
% NE_ridist = nedata.shuffled_stats.NE_reliability_idx;
% NE_ri = mean(nedata.STA_stats.NE_reliability_idx, 2, 'omitnan');
% NE_ridist = [NE_ridist NE_ri];
% NE_ri_pval = sum(NE_ri <= NE_ridist, 2, 'omitnan') ./ size(NE_ridist, 2);
% NE_ri_sig = NE_ri_pval < relidxthresh;
% NE_ri_zscore = nedata.sig_classifying_stats.NE_relidx_zscore;
% NE_ri_sigz = NE_ri_zscore >= zscorethresh;
% 
% neuron_ridist = nedata.shuffled_stats.neuron_reliability_idx;
% neuron_ri = mean(nedata.STA_stats.neuron_reliability_idx, 2, 'omitnan');
% neuron_ridist = [neuron_ridist neuron_ri];
% neuron_ri_pval = sum(neuron_ri <= neuron_ridist, 2, 'omitnan') ./ size(neuron_ridist, 2);
% neuron_ri_sig = neuron_ri_pval < relidxthresh;
% neuron_ri_zscore = nedata.sig_classifying_stats.neuron_relidx_zscore;
% neuron_ri_sigz = neuron_ri_zscore >= zscorethresh;

neuronidx = contains(fieldnames(sig), 'neuron');
temp = struct2cell(sig);
tempz = struct2cell(sigz);

if zscorecheck
    neuron_sig = sum(cell2mat(temp(neuronidx)'), 2) == 3 |...
        sum(cell2mat(tempz(neuronidx)'), 2) >= 2;
    NE_sig = sum(cell2mat(temp(~neuronidx)'), 2) == 3 | ...
        sum(cell2mat(tempz(~neuronidx)'), 2) >= 2;
else
    neuron_sig = sum(cell2mat(temp(neuronidx)'), 2) == 3;
    NE_sig = sum(cell2mat(temp(~neuronidx)'), 2) == 3;
end

% details.NE_ptd_pval = NE_ptd_pval;
% details.NE_moransI_pval = NE_mi_pval;
% details.NE_relidx_pval = NE_ri_pval;
% 
% details.neuron_ptd_pval = neuron_ptd_pval;
% details.neuron_moransI_pval = neuron_mi_pval;
% details.neuron_relidx_pval = neuron_ri_pval;
% 
% details.NE_ptd_zscore = NE_ptd_zscore;
% details.NE_moransI_zscore = NE_mi_zscore;
% details.NE_relidx_zscore = NE_ri_zscore;
% 
% details.neuron_ptd_zscore = neuron_ptd_zscore;
% details.neuron_moransI_zscore = neuron_mi_zscore;
% details.neuron_relidx_zscore = neuron_ri_zscore;

if ~isempty(plotopt)
    
    switch plotopt
        case 'sig'   
            neuron_sta = nedata.sig_stamat;
            NE_sta = nedata.sig_NE_stamat;
        case 'raw'
            neuron_sta = nedata.stamat;
            NE_sta = nedata.NE_stamat;
    end
            
    nf = nedata.nf;
    nlags = nedata.nlags;
    
    cmap = flipud(brewermap(1000, 'RdBu'));
    
    sig_neuron_sta = neuron_sta(neuron_sig,:);
    neuron_num = find(neuron_sig);
    figure;
    colormap(cmap);
    c = 1;
    
    for i = 1:size(sig_neuron_sta, 1)
        
        subplot(4,4,c)
        imagesc(reshape(sig_neuron_sta(i,:), nf, nlags))
        
        cmax = max(sig_neuron_sta(i,:));
        cmin = min(sig_neuron_sta(i,:));
        climit = max(abs([cmax cmin]));
        set(gca, 'clim',[-climit climit])
        title(sprintf('Neuron #%d', neuron_num(i)));
        
        if c == 16 && i ~= size(sig_neuron_sta, 1)
            suptitle('Significant neuronal STAs')
            figure;
            colormap(cmap);
            c = 1;
        else
            c = c+1;
        end
    end
    suptitle('Significant neuronal STAs')

    nsig_neuron_sta = neuron_sta(~neuron_sig, :);
    neuron_num = find(~neuron_sig);
    figure;
    colormap(cmap);
    c = 1;
    
    for i = 1:size(nsig_neuron_sta, 1)
        
        subplot(4,4,c)
        imagesc(reshape(nsig_neuron_sta(i,:), nf, nlags))
        
        cmax = max(nsig_neuron_sta(i,:));
        cmin = min(nsig_neuron_sta(i,:));
        climit = max(abs([cmax cmin]));
        set(gca, 'clim',[-climit climit])
        title(sprintf('Neuron #%d', neuron_num(i)));
        
        if c == 16 && i ~= size(nsig_neuron_sta, 1)
            suptitle('Not significant neuronal STAs')
            figure;
            colormap(cmap);
            c = 1;
        else
            c = c+1;
        end
    end
    suptitle('Not significant neuronal STAs')
    
    sig_NE_sta = NE_sta(NE_sig,:);
    NE_num = find(NE_sig);
    figure;
    colormap(cmap);
    c = 1;
    
    for i = 1:size(sig_NE_sta, 1)
        
        subplot(4,4,c)
        imagesc(reshape(sig_NE_sta(i,:), nf, nlags))
        
        cmax = max(sig_NE_sta(i,:));
        cmin = min(sig_NE_sta(i,:));
        climit = max(abs([cmax cmin]));
        set(gca, 'clim',[-climit climit])
        title(sprintf('NE #%d', NE_num(i)));
        
        if c == 16 && i ~= size(sig_NE_sta, 1)
            suptitle('Significant cNE STAs')
            figure;
            colormap(cmap);
            c = 1;
        else
            c = c+1;
        end
    end
    suptitle('Significant cNE STAs')

    
    nsig_NE_sta = NE_sta(~NE_sig,:);
    NE_num = find(~NE_sig);
    figure;
    colormap(cmap);
    c = 1;
    
    for i = 1:size(nsig_NE_sta, 1) 
        
        subplot(4,4,c)
        imagesc(reshape(nsig_NE_sta(i,:), nf, nlags))
        
        cmax = max(nsig_NE_sta(i,:));
        cmin = min(nsig_NE_sta(i,:));
        climit = max(abs([cmax cmin]));
        set(gca, 'clim',[-climit climit])
        title(sprintf('NE #%d', NE_num(i)));

        if c == 16 && i ~= size(nsig_NE_sta, 1)
            suptitle('Not significant cNE STAs')
            figure;
            colormap(cmap);
            c = 1;
        else
            c = c+1;
        end
    end
    suptitle('Not significant cNE STAs')

end
end
    
    
            
    
    




