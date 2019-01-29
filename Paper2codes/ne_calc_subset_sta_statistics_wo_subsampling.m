function statstruct = ne_calc_subset_sta_statistics_wo_subsampling(exp_site_nedata, spksubset)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);

fn = fieldnames(spksubset);
fn = fn(~cellfun('isempty', regexp(fn, '^spktrain')));

spktrain = cellfun(@(x) {spksubset.(x)}, fn, 'UniformOutput', 0);
spktrain = cell2mat([spktrain{:}]');
stastats = ne_calc_sta_stats_zscore_from_spktrain(spktrain, stimstr.stimulus);

fnprefix = regexp(fn, '(?<=^spktrain_)\S+$', 'match', 'once');
stastatsfn = fieldnames(stastats);

statstruct(length(spksubset)).exp = [];

for i = 1:length(spksubset)
    
    statstruct(i).exp = exp_site_nedata.exp;
    statstruct(i).site = exp_site_nedata.site;
    statstruct(i).stim = exp_site_nedata.stim;
    statstruct(i).sitmlength = exp_site_nedata.stimlength;
    statstruct(i).df = exp_site_nedata.df;
    statstruct(i).fs = exp_site_nedata.fs;
    statstruct(i).depth = exp_site_nedata.depth;
    statstruct(i).probetype = exp_site_nedata.probetype;
    statstruct(i).neuron = spksubset(i).neuron;
    statstruct(i).NE = spksubset(i).ensemble;
    
    for j = 1:length(stastatsfn)
        
        for k = 1:length(fnprefix)
            statstruct(i).([fnprefix{k} '_' stastatsfn{j}]) = stastats.(stastatsfn{j})((k-1)*length(statstruct) + i, :);
        end
    end
end



