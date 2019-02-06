function pvals = ne_calc_silhouette_pvals_spike_subsets(waveforms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

niter = 100;
pvals = zeros(length(waveforms.category), 1);

for i = 1:length(waveforms.category)
    
    fprintf('\nCalculating pvals for unit %d of %d...', i, length(waveforms.category))
    
    cat = waveforms.category{i};
    proj = waveforms.projections{i};
    randsil = zeros(niter, 1);
    
    for j = 1:niter
        
        randcat = cat(randperm(length(cat)));
        randsil(j) = mean(silhouette(proj, randcat));
        
    end
    
    realsil = mean(waveforms.silhouette{i});        
    pvals(i) = sum(realsil <= [realsil; randsil]) ./ (niter+1);

end

fprintf('\n');

