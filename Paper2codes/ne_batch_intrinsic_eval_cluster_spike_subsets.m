function eval_pvals = ne_batch_intrinsic_eval_cluster_spike_subsets(nefiles, evaltype, measuretype)

eval_pvals = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'subset_waveforms')
    if ~exist('subset_waveforms','var')
        continue
    end
    
    fprintf('\nProcessing %s indices for %s...\n', evaltype, nefiles{i})
    eval_pvals{i} = ne_intrinsic_eval_cluster_spike_subsets(subset_waveforms, evaltype, measuretype);
    clear('subset_waveforms')
    
    
end 
    
    