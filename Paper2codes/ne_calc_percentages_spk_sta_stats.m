function stats = ne_calc_percentages_spk_sta_stats(exp_site_nedata, percentages, niter)

% As of 1/9/19, this function takes way too long to run on a single site.
% In any case, one run of this has convinced me thoroughly that spike
% counts make a difference to z-score statistics.

if ~exist('percentages','var')
    percentages = 0.1:0.1:0.9;
end

if ~exist('niter', 'var')
    niter = 20;
end

spktrain = exp_site_nedata.nedata.sta_spktrain;
spkcount = sum(spktrain, 2);
nlags = exp_site_nedata.nedata.nlags;
stimlength = exp_site_nedata.stimlength;
stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);

% ptd = cell(length(percentages), niter);
% moransI = cell(length(percentages), niter);
% RI = cell(length(percentages), niter);

nreps = 100;
pval = 95;

stats(niter).temp = [];

for i = 1:length(percentages)
    
    clc;
    fprintf('Processing %d of %d percentages...\n', i, length(percentages))
    
    for j = 1:niter
        
        tempspk = zeros(size(spktrain));
        
        for k = 1:length(spkcount)
            
            tempspk(k,:) = sub_sample_spktrain(spktrain(k,:), ...
                round(spkcount(k) - percentages(i) * spkcount(k)));
            
        end
        
        [realval.ptd, realval.moransI, realval.RI, sta] = ne_calc_STA_statistics_from_spktrain(...
            tempspk, stimstr.stimulus, nlags, stimlength);
        
        
        [~, shuffval.ptd, shuffval.moransI, shuffval.RI] =...
            ne_calc_sig_sta_statistics (sta, tempspk, stimstr.stimulus,...
            nreps, nlags, stimlength, pval);
        
        fn = fieldnames(realval);
        
        for ii = 1:length(fn)
            
            % Get z-score
            tempS = shuffval.(fn{ii});
            mu = nanmean(tempS, 2);
            nansum = sum(isnan(tempS), 2);
            sigma = nanstd(tempS, 0, 2);
            
            tempR = realval.(fn{ii});
            tempR = nanmean(tempR, 2);
            
            stats(j).([fn{ii} '_zscore_' num2str(percentages(i)*100)]) = (tempR - mu) ./ sigma;
            
            % Get p-val
            dist = [tempS tempR];
            stats(j).([fn{ii} '_pval_' num2str(percentages(i)*100)]) = (sum(dist...
                >= tempR, 2) - nansum) ./ (size(tempS, 2) - nansum);
        end                                    
    end
end
        
stats = rmfield(stats, 'temp');        
              
    
end

