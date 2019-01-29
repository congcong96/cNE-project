function sig_sta = ne_group_classified_STAs_based_on_NE(exp_site_nedata)
%UNTITLED23 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;

NEmembers = nedata.NEmembers;
sig_neuron = nedata.sig_neuron_sta;
sig_NE = nedata.sig_NE_sta;

sig_sta(length(NEmembers)).NE = [];

for i = 1:length(NEmembers)
    
    sig_sta(i).NE = i;
    sig_sta(i).NEmembers = NEmembers{i};
    sig_sta(i).NE_sig = sig_NE(i);
    sig_sta(i).neuron_sig = sig_neuron(NEmembers{i});

end

