function numspk = ne_batch_find_number_of_spikes_per_NE_event(nefiles, NEthreshalpha, stat)

numspk = cell(length(nefiles), 1);

for j = 1:length(nefiles)
    
    load(nefiles{j}, 'exp_site_nedata')

    nedata = exp_site_nedata.nedata;
    spktrain = nedata.spktrain;
    NEact = nedata.Activities;

    NEtidx = nedata.NEthresh_alpha == NEthreshalpha; 
    NEthresh = nedata.NEthresh(NEtidx,:);
    NEmembers = nedata.NEmembers;
    NEsize = cellfun('length', NEmembers);

    numspk{j} = zeros(length(NEthresh),1);

    for i = 1:length(NEthresh)

        idx = (NEact(i,:) >= NEthresh(i));
        tempspktrain = logical(spktrain(NEmembers{i}, idx));
        switch stat
            case 'min'
                numspk{j}(i) = min(sum(tempspktrain,1));
            case 'mean'
                numspk{j}(i) = mean(sum(tempspktrain,1));
            case 'max'
                numspk{j}(i) = max(sum(tempspktrain,1));
            case 'mode'
                numspk{j}(i) = mode(sum(tempspktrain,1));
            case 'median'
                numspk{j}(i) = median(sum(tempspktrain,1));
            case 'minperc'
                numspk{j}(i) = min(sum(tempspktrain,1)) / NEsize(i);
            case 'meanperc'
                numspk{j}(i) = mean(sum(tempspktrain,1)) / NEsize(i);
            case 'maxperc'
                numspk{j}(i) = max(sum(tempspktrain,1)) / NEsize(i);
            case 'modeperc'
                numspk{j}(i) = mode(sum(tempspktrain,1)) / NEsize(i);
            case 'medianperc'
                numspk{j}(i) = median(sum(tempspktrain,1)) / NEsize(i);
            otherwise
                error('Invalid stat option')
        end

    end
    
end

