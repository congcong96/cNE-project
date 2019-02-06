function pvals = ne_calc_F1_pvals_spike_subsets(waveforms)

kiter = 10;
niter = 100;
pvals = zeros(length(waveforms.category), 1);

for i = 1:length(waveforms.category)
    
    fprintf('\nCalculating pvals for unit %d of %d...', i, length(waveforms.category))
    
    cat = waveforms.category{i};
    proj = waveforms.projections{i};
    
    cat = cat(:);
    
    uniqcat = unique(cat);
    
    if ismember(0, uniqcat)
        cat = cat + 1;
        uniqcat = uniqcat + 1;
    end
    
    randF = zeros(niter, kiter);
    realF = zeros(1, kiter);
    
    for k = 1:kiter
        randpred = kmeans(proj, length(uniqcat), 'MaxIter', 10000);

        for j = 1:niter
            randcat = cat(randperm(length(cat)));
            randF(j, k) = mean(calc_F_for_clusters(randcat, randpred, 1));
        end

        realF(k) = mean(calc_F_for_clusters(cat, randpred, 1));        
        
    end
    
    meanF = mean(realF);
    pvals(i) = sum(meanF <= [meanF; randF(:)]) ./ (niter*kiter+1);
    
    if pvals(i) == 0
        keyboard
    end
    
end

fprintf('\n');

end

