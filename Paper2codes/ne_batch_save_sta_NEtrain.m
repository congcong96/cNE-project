function ne_batch_save_sta_NEtrain(nefiles, NEthreshalpha)

if ~exist('NEthreshalpha','var')
    NEthreshalpha = 99.5;
end

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    NEact = nedata.Activities;
    NEmembers = nedata.NEmembers;
    spktrain = nedata.sta_spktrain;
    
    NEtidx = nedata.NEthresh_alpha == NEthreshalpha;
    
    NEthresh = nedata.NEthresh(NEtidx,:);
    
    exp_site_nedata.nedata.sta_NEtrain = ne_upsample_NEact_using_member_neuron_activity(NEact, NEmembers, spktrain, NEthresh);
    
    save(nefiles{i}, 'exp_site_nedata', '-append');

end
    