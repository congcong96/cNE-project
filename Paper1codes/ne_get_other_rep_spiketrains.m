function othertrains = ne_get_other_rep_spiketrains(files)


for i = 1:length(files)
    load(files{i})
    othertrains{i} = exp_site_nedata.nedata.spktrain;

end

