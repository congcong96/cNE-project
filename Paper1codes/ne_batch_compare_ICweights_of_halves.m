function [corrmat, IC1, IC2] = ne_batch_compare_ICweights_of_halves(files)

for i = 1:length(files)
    load(files{i}, 'spktrain')
    
    [corrmat{i}, IC1{i}, IC2{i}] = ne_compare_ICweights_of_halves(spktrain, 20, 'interleaved');
    
    clear('spktrain')
end

