function corrvec = ne_get_STAsim_of_member_STAs(sig_sta, sigopt)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

corrvec = cell(length(sig_sta), 1);

for i = 1:length(sig_sta)
    
    if sigopt
        neuron_sta = sig_sta(i).neuron_sig_sta(sig_sta(i).sig_neurons,:);
    else
        neuron_sta = [sig_sta(i).neuron_sta];
    end
    if size(neuron_sta, 1) == 1
        continue
    else
        sta_corr = corr(neuron_sta');
    end
    
    corrvec{i} = sta_corr(triu(sta_corr, 1) ~= 0);    

end