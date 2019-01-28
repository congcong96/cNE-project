function [rate] = ne_calc_neuron_firing_rate(exp_site_nedata)


df = exp_site_nedata.df;
nedata = exp_site_nedata.nedata;
spktrain = logical(nedata.spktrain);
len = size(spktrain, 2);

tottime = 0.5 * df * len / 1000; % in seconds

rate = sum(spktrain,2)./tottime;