function sig_sta_subset = ne_calc_sig_sta_subset_struct(sig_sta, calc_neurons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

numneurons = sum(cell2mat({sig_sta.sig_neurons}') == calc_neurons);
sig_sta_subset(numneurons).filename = [];
c = 1;

for i = 1:length(sig_sta)
    
    load(sig_sta(i).filename, 'exp_site_nedata')
    NE = sig_sta(i).NE;
    nedata = exp_site_nedata.nedata;
    NEtrain = nedata.sta_NEtrain(NE,:);
    
    members = nedata.NEmembers{NE};
    spktrain = nedata.spktrain(members, :);

    NEsubset = ne_get_NEsubsets_from_member_spiketrain(NEtrain, spktrain, 'individual_members');        
    
    for j = 1:length(sig_sta(i).neurons)        
              
        if sig_sta(i).sig_neurons(j) == calc_neurons
           
            sig_sta_subset(c).filename = sig_sta(i).filename;
            sig_sta_subset(c).NE = sig_sta(i).NE;
            sig_sta_subset(c).neuron = sig_sta(i).neurons(j);
            sig_sta_subset(c).all_spktrain = nedata.sta_spktrain(sig_sta(i).neurons(j),:);         
            sig_sta_subset(c).wNE_subset_spktrain = NEsubset(j,:);
            sig_sta_subset(c).woNE_subset_spktrain = sig_sta_subset(c).all_spktrain - sig_sta_subset(c).wNE_subset_spktrain;
            sig_sta_subset(c).all_count = sum(sig_sta_subset(c).all_spktrain);
            sig_sta_subset(c).wNE_count = sum(sig_sta_subset(c).wNE_subset_spktrain);
            sig_sta_subset(c).woNE_count = sum(sig_sta_subset(c).woNE_subset_spktrain);
            c = c+1;
        else
            continue
        end
    end
end
end

