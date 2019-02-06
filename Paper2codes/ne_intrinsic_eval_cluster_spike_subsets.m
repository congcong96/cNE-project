function pvals = ne_intrinsic_eval_cluster_spike_subsets(waveforms, evaltype, measuretype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

niter = 1000;
pvals = zeros(length(waveforms.category), 1);

for i = 1:length(waveforms.category)
    
    fprintf('\nCalculating pvals for unit %d of %d...', i, length(waveforms.category))
    
    cat = waveforms.category{i};
    proj = waveforms.projections{i};
    randeval = zeros(niter, 1);
    
    for j = 1:niter
        
        randcat = cat(randperm(length(cat)));
        switch evaltype
            case 'Dunn'
                randeval(j) = mean(calc_Dunn_like_index(proj, randcat, measuretype));                
            case 'DaviesBouldin'
                randeval(j) = calc_DaviesBouldin_index(proj, randcat, measuretype);
        end         
        
    end
    
    switch evaltype
        case 'Dunn'
            realeval = mean(calc_Dunn_like_index(proj, cat, measuretype));
            pvals(i) = sum(realeval <= [realeval; randeval]) ./ (niter+1);
        case 'DaviesBouldin'
            realeval = calc_DaviesBouldin_index(proj, cat, measuretype);
            pvals(i) = sum(realeval >= [realeval; randeval]) ./ (niter+1);
    end              
    

end

fprintf('\n');

