function ne_batch_calc_subsets_repeats_NE_stats(files, niter, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


for i = 1:niter
    
    if exist(sprintf('rep2NE%d.mat',i), 'file')
        fprintf('\nrep2NE%d.mat already exists! Skipping... \n', i)
        continue
    end

    clc;
    fprintf('Iteration %d of %d\n',i,niter)
    if isempty(varargin)
        repNE = ne_calc_subsets_repeats_NE_stats(files);
    else
        repNE = ne_calc_subsets_repeats_NE_stats(files, 'goodneurons', varargin{1});
    end
    
    save(sprintf('repNE%d',i), 'repNE')    
    clear('repNE')
    

end

