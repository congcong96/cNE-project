function ne_batch_classify_significant_STAs(nefiles)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

pvalthresh.ptd_thresh = 1;
pvalthresh.moransI_thresh = 1;
pvalthresh.reliability_idx_thresh = 0.05;

for i = 1:length(nefiles)
    
    fprintf('\nClassifying significant STAs for %s...\n', nefiles{i})
    
    load(nefiles{i}, 'exp_site_nedata')

    [neuron_sig, NE_sig, details] = ne_classify_significant_STAs(exp_site_nedata,...
        'pvalthresh', pvalthresh);
    
    exp_site_nedata.nedata.sig_neuron_sta = neuron_sig;
    exp_site_nedata.nedata.sig_NE_sta = NE_sig;
    exp_site_nedata.nedata.sig_classifying_stats = details;
        
    save(nefiles{i}, 'exp_site_nedata', '-append');
    clear('exp_site_nedata')


end

