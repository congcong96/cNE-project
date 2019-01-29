function [bvec, pval] = ne_plot_MIdiff_vs_ICweight(exp_site_nedata, NEneuroninfo, plotopt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

num_member_thresh = 5;

NEneuroninfo(cellfun('isempty', {NEneuroninfo.fraction})) = [];

if isempty(NEneuroninfo)
    bvec = [];
    pval = [];
    return
end

NEs = [NEneuroninfo.NE];
[count, val] = rude(NEs);

NEs_toplot = val(count >= num_member_thresh); 

if isempty(NEs_toplot)
    bvec = [];
    pval = [];
    return
end

pval = zeros(length(NEs_toplot), 1);
bvec = zeros(length(NEs_toplot), 2);

for i = 1:length(NEs_toplot)
    
    tempNEneuron = NEneuroninfo([NEneuroninfo.NE] == NEs_toplot(i));
    infodiff = [tempNEneuron.neuron_info_extrap] - [tempNEneuron.NE_info_extrap];
    
    neurons = [tempNEneuron.neuron];
    ICwt = exp_site_nedata.nedata.Patterns(neurons, NEs_toplot(i));
    
    if sum(ICwt<0) > sum(ICwt>0)
        ICwt = -ICwt;
    end
    
    [b,~,~,~,stats] = regress(infodiff(:), [ones(length(ICwt), 1), ICwt(:)]);
    
    pval(i) = stats(3);
    bvec(i, :) = b';
    
    if plotopt      
        figure;
        scatter(ICwt, infodiff); 
        h = refline(b(2), b(1));
        h.Color = 'r';
    end


end

