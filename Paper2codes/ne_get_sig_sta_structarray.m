function sig_sta = ne_get_sig_sta_structarray(nefiles)

sig_sta(length(nefiles) * 10).filename = [];
c = 1;

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    NEmembers = nedata.NEmembers;
        
    for j = 1:length(NEmembers)
        
        sig_sta(c).filename = nefiles{i};
        sig_sta(c).IC_weights = nedata.Patterns(:,j);
        sig_sta(c).CI = nedata.CI;
        sig_sta(c).NE = j;
        sig_sta(c).neurons = NEmembers{j};
        sig_sta(c).num_members = length(NEmembers{j});
        sig_sta(c).sig_NE = nedata.sig_NE_sta(j);
        sig_sta(c).sig_neurons = nedata.sig_neuron_sta(NEmembers{j});
        sig_sta(c).NE_event_count = sum(nedata.sta_NEtrain(j,:));
        sig_sta(c).neuron_spike_count = sum(nedata.sta_spktrain(NEmembers{j},:), 2);
        sig_sta(c).NE_sta = nedata.NE_stamat(j,:);
        sig_sta(c).neuron_sta = nedata.stamat(NEmembers{j},:);
        sig_sta(c).NE_sig_sta = nedata.sig_NE_stamat(j,:);
        sig_sta(c).neuron_sig_sta = nedata.sig_stamat(NEmembers{j}, :);
        sig_sta(c).percentage_sig_neurons = sum(sig_sta(c).sig_neurons)/length(sig_sta(c).sig_neurons);
        sig_sta(c).NE_ptd_pval = nedata.sig_classifying_stats.NE_ptd_pval(j);
        sig_sta(c).neuron_ptd_pval = nedata.sig_classifying_stats.neuron_ptd_pval(NEmembers{j});
        sig_sta(c).NE_moransI_pval = nedata.sig_classifying_stats.NE_moransI_pval(j);
        sig_sta(c).neuron_moransI_pval = nedata.sig_classifying_stats.neuron_moransI_pval(NEmembers{j});
        sig_sta(c).NE_reliability_idx_pval = nedata.sig_classifying_stats.NE_reliability_idx_pval(j);
        sig_sta(c).neuron_reliability_idx_pval = nedata.sig_classifying_stats.neuron_reliability_idx_pval(NEmembers{j});
        sig_sta(c).NE_ptd_zscore = nedata.sig_classifying_stats.NE_ptd_zscore(j);
        sig_sta(c).neuron_ptd_zscore = nedata.sig_classifying_stats.neuron_ptd_zscore(NEmembers{j});
        sig_sta(c).NE_moransI_zscore = nedata.sig_classifying_stats.NE_moransI_zscore(j);
        sig_sta(c).neuron_moransI_zscore = nedata.sig_classifying_stats.neuron_moransI_zscore(NEmembers{j});
        sig_sta(c).NE_reliability_idx_zscore = nedata.sig_classifying_stats.NE_reliability_idx_zscore(j);
        sig_sta(c).neuron_reliability_idx_zscore = nedata.sig_classifying_stats.neuron_reliability_idx_zscore(NEmembers{j});
        
        if sig_sta(c).sig_NE == 0 && sig_sta(c).percentage_sig_neurons > 0.5
            sig_sta(c).category = 'destructive';
        elseif sig_sta(c).sig_NE == 0 && sig_sta(c).percentage_sig_neurons <= 0.5
            sig_sta(c).category = 'stimindependent';
        elseif sig_sta(c).sig_NE == 1 && sig_sta(c).percentage_sig_neurons < 0.5
            sig_sta(c).category = 'constructive';
        elseif sig_sta(c).sig_NE == 1 && sig_sta(c).percentage_sig_neurons >= 0.5
            sig_sta(c).category = 'facilitative';
        else
            error('Something''s wrong!')
        end
        
        c = c+1; 
    end
    
    sig_sta(cellfun('isempty', {sig_sta.NE})) = [];
    
end

sig_sta([sig_sta.num_members] == 1) = [];
   
