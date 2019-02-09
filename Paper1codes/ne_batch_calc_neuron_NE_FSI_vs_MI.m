function compcell = ne_batch_calc_neuron_NE_FSI_vs_MI(files)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
 
[FSI_ne, FSI_neuron] = ne_batch_calc_neuron_NE_FSI(files);

compcell = cell(length(files), 1);

for i = 1:length(files)

    load(files{i}, 'NEneuroninfo2')
    compstruct(length(NEneuroninfo2)).neuron = [];

    for j = 1:length(compstruct)
        
        compstruct(j).neuron = NEneuroninfo2(j).neuron;
        compstruct(j).NE = NEneuroninfo2(j).NE;
        
        if ~isempty(NEneuroninfo2(j).NE_info_extrap)
            compstruct(j).MIdiff = NEneuroninfo2(j).NE_info_extrap - NEneuroninfo2(j).neuron_info_extrap;
            compstruct(j).FSIdiff = FSI_ne{i}(compstruct(j).NE) - FSI_neuron{i}(compstruct(j).neuron);  
        end
        
    end
    
    compcell{i} = compstruct;
    clear('compstruct')
       

end




