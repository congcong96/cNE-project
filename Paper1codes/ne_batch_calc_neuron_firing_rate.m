function rate = ne_batch_calc_neuron_firing_rate(NEfiles)

rate = cell(1, length(NEfiles));

for i = 1:length(NEfiles)
    load(NEfiles{i}, 'exp_site_nedata')
    rate{i} = ne_calc_neuron_firing_rate(exp_site_nedata);
end

end

