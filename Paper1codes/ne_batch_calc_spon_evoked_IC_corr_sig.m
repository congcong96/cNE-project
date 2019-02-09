function corrvalstruct = ne_batch_calc_spon_evoked_IC_corr_sig(corrvalstruct)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


for i = 1:length(corrvalstruct)
    
    absshuffcorr = abs(corrvalstruct(i).shuffledcorrvals);
    
    sponpvals = zeros(1,length(corrvalstruct(i).sponcorrvals));
    evokedpvals = zeros(1,length(corrvalstruct(i).evokedcorrvals));
    
    for j = 1:length(corrvalstruct(i).sponcorrvals)
        
        sponpvals(j) = (sum(abs(corrvalstruct(i).sponcorrvals(j)) <= absshuffcorr)...
            + 1) ./ (length(absshuffcorr) + 1);
        
    end
    
    for k = 1:length(corrvalstruct(i).evokedcorrvals)
        

        
        evokedpvals(k) = (sum(abs(corrvalstruct(i).evokedcorrvals(k)) <= ...
            absshuffcorr) + 1) ./ (length(absshuffcorr) + 1);
    end
    
    corrvalstruct(i).sponpvals = sponpvals;
    corrvalstruct(i).evokedpvals = evokedpvals;       
                
        

end

