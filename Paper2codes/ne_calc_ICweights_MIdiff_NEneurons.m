function ICwtMIstruct = ne_calc_ICweights_MIdiff_NEneurons(exp_site_nedata, NEneuroninfo)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NEneuroninfo(cellfun('isempty', {NEneuroninfo.fraction})) = [];

if isempty(NEneuroninfo)
    ICwtMIstruct = [];
    return
end

NEs = [NEneuroninfo.NE];
uniqueNEs = unique(NEs);

ICwtMIstruct(length(uniqueNEs)).NE = [];

for j = 1:length(uniqueNEs)
    
    tempNEneuron = NEneuroninfo([NEneuroninfo.NE] == uniqueNEs(j));
    
    if length(tempNEneuron) == 1
        ICwtMIstruct(j).NE = [];
        ICwtMIstruct(j).neurons = [];
        ICwtMIstruct(j).sig_NE = [];
        ICwtMIstruct(j).sig_neurons = [];
        ICwtMIstruct(j).MIdiff = [];
        ICwtMIstruct(j).ICwt = [];
        continue
    end    
    
    infodiff = [tempNEneuron.neuron_info_extrap] - [tempNEneuron.NE_info_extrap];
    
    neurons = [tempNEneuron.neuron];
    ICwt = exp_site_nedata.nedata.Patterns(neurons, uniqueNEs(j));
    
    if sum(ICwt<0) > sum(ICwt>0)
        ICwt = -ICwt;
    elseif sum(ICwt<0) == sum(ICwt>0)
        [~, idx] = max(abs(infodiff));
        if ICwt(idx) < 0
            ICwt = -ICwt;
        end        
    end
    
%     infodiff = infodiff(ICwt > 0);
    
%     if length(infodiff) <= 1
%         continue
%     end
    
    ICwtMIstruct(j).NE = uniqueNEs(j);
    ICwtMIstruct(j).neurons = neurons;
    ICwtMIstruct(j).sig_NE = exp_site_nedata.nedata.sig_NE_sta(uniqueNEs(j));
    ICwtMIstruct(j).sig_neurons = exp_site_nedata.nedata.sig_neuron_sta(neurons);
    ICwtMIstruct(j).MIdiff = infodiff;
    ICwtMIstruct(j).ICwt = ICwt;

    
%     ICwt = ICwt(ICwt > 0);
    
%     [~, maxidx] = max(infodiff);
%     [~, minidx] = min(infodiff);
    
  
%     ICwtMIstruct(j).minMI_ICwt = ICwt(minidx);
    
%     [~, maxidx] = max(ICwt);
%     [~, minidx] = min(ICwt);
%      
%     ICwtMIstruct(j).maxICwt_MI = infodiff(maxidx);
%     ICwtMIstruct(j).minICwt_MI = infodiff(minidx);    
%     
%     [~, maxidx] = max(diffinfo(tempidx));
%     [~, minidx] = min(diffinfo(tempidx));
%     maxneuron = NEneuroninfo(tempidx(maxidx)).neuron;
%     minneuron = NEneuroninfo(tempidx(minidx)).neuron;
%     tempmax = exp_site_nedata.nedata.Patterns(maxneuron, uniqueNEs(j));
%     tempmin = exp_site_nedata.nedata.Patterns(minneuron, uniqueNEs(j));
%     
%     if tempmax < 0 && tempmin < 0
%         tempmax = abs(tempmax);
%         tempmin = abs(tempmin);
%     elseif tempmax < 0 && tempmin > 0
%         tempmax = -tempmax;
%         tempmin = -tempmin;
%     end
%     
%     max_ICwt(j) = tempmax;
%     min_ICwt(j) = tempmin;
    
end
    
    

