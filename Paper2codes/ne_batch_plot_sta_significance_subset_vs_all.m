function [p, varargout] = ne_batch_plot_sta_significance_subset_vs_all(sig_sta, paramopt, spkthresh)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('paramopt', 'var')
    paramopt = 1;
end

if ~exist('spkthresh', 'var')
    spkthresh = 0;
end

sig_all = zeros(length(sig_sta), 1);
sig_subset = zeros(length(sig_sta), 1);

subset_num_sig = cell(length(sig_sta), 1);

for i = 1:length(sig_sta)
    
    load(sig_sta(i).filename, 'nosample_subsetstats');
    
    NE = sig_sta(i).NE;
    neurons = sig_sta(i).neurons;
    
    neuronidx = logical(sum(cell2mat(arrayfun(@(x) [nosample_subsetstats.neuron] == x, neurons, 'UniformOutput' ,0)), 1));
    NEidx = [nosample_subsetstats.NE] == NE;
    wNE_count = [nosample_subsetstats.w_NE_count];
    
    idx = neuronidx & NEidx & wNE_count > spkthresh;
    
%     if sum(idx) == 0
%         continue
%     end
    
    wNE_ri_pval = [nosample_subsetstats(idx).w_NE_reliability_idx_pval];
    sig_subset(i) = sum(wNE_ri_pval < 0.05) / length(wNE_ri_pval);
    all_ri_pval = [nosample_subsetstats(idx).all_reliability_idx_pval];
    sig_all(i) = sum(all_ri_pval < 0.05) / length(all_ri_pval);

    subset_num_sig{i} = wNE_ri_pval < 0.05;
    
end
    
p = plot_paired_data([sig_all sig_subset], {'All spikes', 'cNE subsets'}, 1, paramopt);
ylabel('Ratio of significant STRFs')
print_mfilename(mfilename);

varargout{1} = sig_all;
varargout{2} = sig_subset;
varargout{3} = subset_num_sig;

