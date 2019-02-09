function surrnedata = ne_shuffled_rnrep_trials_NE_data_processing(spk, trigger, binsize, niter)

surrnedata(niter).orderidx = [];

for i = 1:niter
    [shuffledspkmatrix, orderidx] = ne_shuffle_rnrep_trial_order(spk, trigger, binsize);
    
    nedata = ne_data_processing_for_shuffled_NE_raster(shuffledspkmatrix, binsize);
    surrnedata(i).nedata = nedata;
    surrnedata(i).orderidx = orderidx;
    close all
    
end
