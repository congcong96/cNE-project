function [spkcount, NEeventcount] = ne_batch_calc_num_spikes_and_NEevents(nefiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load(nefiles{1})

alpha = exp_site_nedata.nedata.NEthresh_alpha;

spkcount = cell(length(nefiles),1);
NEeventcount = cell(length(nefiles), length(alpha));

for i = 1:length(nefiles)
    
    load(nefiles{i})
    nedata = exp_site_nedata.nedata;
    spkcount{i} = sum(nedata.spktrain, 2);
    NEact = nedata.Activities;
    
    for j = 1:length(alpha)         
        nethresh = exp_site_nedata.nedata.NEthresh(j,:);
        NEraster = zeros(size(NEact));
        for k = 1:length(nethresh)
            NEraster(k,:) = NEact(k,:) >= nethresh(k);           
        end
        NEeventcount{i,j} = sum(NEraster,2);
    end          

end

