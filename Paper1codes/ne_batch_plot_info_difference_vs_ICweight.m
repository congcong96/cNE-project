function [ICwt, infodiff] = ne_batch_plot_info_difference_vs_ICweight(nefiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

infodiff = cell(length(nefiles),1);
ICwt = cell(length(nefiles),1);

for j = 1:length(nefiles)
    
    temp = load(nefiles{j}, 'NEneuroninfo2', 'exp_site_nedata');
    NEneuroninfo = temp.NEneuroninfo2;
    exp_site_nedata = temp.exp_site_nedata;

    infodiff{j} = zeros(length(NEneuroninfo),1);
    ICwt{j} = zeros(length(NEneuroninfo),1);

    for i = 1:length(NEneuroninfo)
        
        if isempty(NEneuroninfo(i).NE_info_extrap)
            continue
        end

        neuron = NEneuroninfo(i).neuron;
        NE = NEneuroninfo(i).NE;

        infodiff{j}(i) = NEneuroninfo(i).NE_info_extrap - NEneuroninfo(i).neuron_info_extrap;
        ICwt{j}(i) = abs(exp_site_nedata.nedata.Patterns(neuron, NE));

    end

end

infodiff = cell2mat(infodiff);
ICwt = cell2mat(ICwt);

zeroidx = ICwt == 0;
ICwt(zeroidx) = [];
infodiff(zeroidx) = [];

figure;

scatter(ICwt, infodiff, 10, 'filled');
xlabel('IC weight')
ylabel('Info NE - info neuron (bits/spike)')

keyboard()