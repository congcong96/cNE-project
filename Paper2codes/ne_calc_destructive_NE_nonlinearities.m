function fiofit = ne_calc_destructive_NE_nonlinearities(des_sigsta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(des_sigsta)
    
    load(des_sigsta(i).filename, 'exp_site_nedata')
    sig_neurons = des_sigsta(i).neurons(des_sigsta(i).sig_neurons);
    cNE = des_sigsta(i).NE;
    
    stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);    
    nedata = exp_site_nedata.nedata;
    
    for j = 1:length(sig_neurons)
        
        [xprior, xposterior] = ne_sta_stimulus_projection(reshape(...
            nedata.stamat(sig_neurons(j),:), nedata.nf, nedata.nlags),...
            nedata.sta_NEtrain(cNE,:), stimstr.stimulus);
        [px,pspk,pxspk,xbinedges] = calc_px_pspk_pxspk(xprior,xposterior);
        pspkx = pspk .* pxspk ./ px;
        xbins = edge2center(xbinedges);    
        fiofit{i}(j) = js_hvv_fio_fit(xbins(:), pspkx, pspk); 
        
    end


end

