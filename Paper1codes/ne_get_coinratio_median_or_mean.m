function ratio = ne_get_coinratio_median_or_mean(allratio, mtype)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
fn = fieldnames(allratio);

for i = 1:length(fn)
    ratiocell = allratio.(fn{i});
    
    if strcmp(fn{i},'repeat')
        switch mtype
            case 'mean'
                temp = cellfun(@mean, ratiocell, 'UniformOutput',0);
                ratio.(fn{i}) = cellfun(@(x) x(:), temp, 'UniformOutput',0);
            case 'median'
                temp = cellfun(@median, ratiocell, 'UniformOutput', 0);
                ratio.(fn{i}) = cellfun(@(x) x(:), temp, 'UniformOutput',0);
        end
    else
        switch mtype
            case 'mean'
                for j = 1:size(ratiocell,1)
                    temp = mean(cell2mat(ratiocell(j,:)), 2);
                    ratio.(fn{i}){j} = temp(:);                     
                end
            case 'median'
                for j = 1:size(ratiocell,1)
                    temp = median(cell2mat(ratiocell(j,:)), 2);
                    ratio.(fn{i}){j} = temp(:);
                end
        end
    end
end



end

