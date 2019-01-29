function [maxidx, bounds] = ne_batch_find_sig_freq_range_for_MID(MID_NEs)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

nefiles = unique({MID_NEs.filename});

c = 1;
maxidx = cell(length(MID_NEs),1);
bounds = cell(length(MID_NEs),1);

for i = 1:length(nefiles)
    
    NEs = [MID_NEs(strcmp(nefiles{i}, {MID_NEs.filename})).NE];
    load(nefiles{i}, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    nf = nedata.nf;
    nlags = nedata.nlags;
    NEmembers = nedata.NEmembers;
    
    for j = 1:length(NEs)
        
        NEstamat = reshape(nedata.NE_stamat(NEs(j),:), nf, nlags);
        [maxidx{c}(1), bounds{c}(1,:)] = find_most_sig_freq_range_for_MID(NEstamat, 3);
        
        for k = 1:length(NEmembers{NEs(j)})
            
            neuronstamat = reshape(nedata.stamat(NEmembers{NEs(j)}(k),:), nf, nlags);
            [maxidx{c}(k+1), bounds{c}(k+1,:)] = find_most_sig_freq_range_for_MID(neuronstamat, 3);
            
        end
        c = c+1;
    end
                      


end

