function neactcount = ne_calc_num_NE_events(exp_site_nedata)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;
neact = nedata.Activities;
nethresh = nedata.NEthresh;

neactbaseline = mode(neact, 2);
neactcount = zeros(size(nethresh));

for i = 1:size(nethresh, 1)
        
    for j = 1:size(nethresh, 2)
        
        if nethresh(i,j) < neactbaseline(j)            
            neactcount(i,j) = inf;            
        else
            neactcount(i,j) = sum(neact(j,:) >= nethresh(i,j));
        end
    end
end
            
        

end

