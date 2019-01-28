function rtfcell = ne_batch_calc_with_or_without_NE_rtf(nefiles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

rtfcell = cell(length(nefiles),1);

for i = 1:length(nefiles)
    
    load(nefiles{i})
    fprintf('\nCalculating RTF values for file %d of %d...',i, length(nefiles))
    rtfcell{i} = ne_plot_with_or_without_assembly_rtf(exp_site_nedata);

end

fprintf('\n')
