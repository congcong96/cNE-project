function surrcoupratio = ne_calc_model_coupling_ratio(files, varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



surrcoupratio = cell(length(files),1);

for i = 1:length(files)
    
    fprintf('\nProcessing file %d of %d...', i, length(files)) 
    load(files{i})
    fn = fieldnames(surrspktrain);
    numneurons = size(surrspktrain.(fn{1}){1}, 1);
    surrcoupratio{i} = zeros(numneurons, length(surrspktrain.(fn{1})));
    
    for j = 1:length(surrspktrain.(fn{1}))
        if ~isempty(varargin)
            surrcoupratio{i}(:,j) = calc_coupling_ratio_with_pfr(1, surrspktrain.(fn{1}){j}, varargin{1});
        else
            surrcoupratio{i}(:,j) = calc_coupling_ratio_with_pfr(1, surrspktrain.(fn{1}){j});
        end

    end


end


fprintf('\n')