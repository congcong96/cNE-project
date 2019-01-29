function ne_batch_calc_num_NE_events(nefiles)

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')  
    necount = ne_calc_num_NE_events(exp_site_nedata);
    exp_site_nedata.nedata.NEeventcount = necount;
    save(nefiles{i}, 'exp_site_nedata', '-append')
    clear('exp_site_nedata')
end