function spktimesub = ne_get_spktime_subset_from_spktrain_subset(spktimes, spktrain, edges)

assert(size(spktrain,2) == length(edges) - 1)

ledge = edges(1:end-1);
redge = edges(2:end);
spktimesub = zeros(1, sum(spktrain));

c = 1;

for i = 1:length(spktrain)
    
    if spktrain(i) > 0
        
        temp = spktimes(spktimes >= ledge(i) & spktimes < redge(i));
        if length(temp) > spktrain(i)
            temp = randsample(temp,spktrain(i));            
        end
        spktimesub(c:c+length(temp)-1) = temp;
        c = c+length(temp);
    end

end

