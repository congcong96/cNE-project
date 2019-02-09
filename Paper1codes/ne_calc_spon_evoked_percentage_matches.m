function [perc, sigcorrvals] = ne_calc_spon_evoked_percentage_matches(corrvalstruct)


folder = 'I:\Cell_Assemblies\MountainSortNEs\sponcomparisons';

% sponmatch = zeros(length(corrvalstruct),1);
% sponlength = zeros(length(corrvalstruct),1);
% evokedmatch = zeros(length(corrvalstruct),1);
% evokedlength = zeros(length(corrvalstruct),1);
sigcorrvals = cell(length(corrvalstruct),1);
sigevoked = cell(length(corrvalstruct),1);
sigspon = cell(length(corrvalstruct),1);

for i = 1:length(corrvalstruct)
    
    filepath = fullfile(folder, corrvalstruct(i).file);
    load(filepath,'exp_site_nedata')
    
    nedata = exp_site_nedata.nedata;
    sponIC = nedata.spon_patterns;
    evokedIC = nedata.evoked_patterns;
    
    corrmat = abs(corr(sponIC,evokedIC));
    
    shuffcorrvals = corrvalstruct(i).shuffledcorrvals;
    thresh = prctile(shuffcorrvals, 99);
    
    %alternative method
    maxevoked = max(corrmat, [], 1);
    maxspon = max(corrmat, [], 2);
    
    sigevoked{i} = maxevoked > thresh;
    sigspon{i} = maxspon > thresh;
    
    
    sigcorrvals{i} = unique([maxevoked(sigevoked{i})'; maxspon(sigspon{i})]);
    
    
%     original method    
 
%     
%     sigidx = corrmat > thresh;
%     
%     sponmatch(i) = sum(sum(sigidx,2) >= 1);
%     sponlength(i) = size(sigidx,1);
%     
%     evokedmatch(i) = sum(sum(sigidx,1) >= 1);
%     evokedlength(i) = size(sigidx,2);
%     
%     sigcorrvals{i} = corrmat(sigidx);

end

% totalNEs = sum(sponlength) + sum(evokedlength);
sigevoked = cell2mat(sigevoked')';
sigspon = cell2mat(sigspon);
totalNEs = length(sigevoked) + length(sigspon);
% totalmatched = sum(sponmatch) + sum(evokedmatch);
totalmatched = sum(sigevoked) + sum(sigspon);
perc = totalmatched/totalNEs;

sigcorrvals = cell2mat(sigcorrvals);
edges = 0:0.05:1;

figure;
histogram(sigcorrvals, edges, 'Normalization','probability');
xlabel('Correlation value')
ylabel('Ratio')
y = ylim;

text(0.1, y(2) * 0.6, sprintf('n = %d', length(sigcorrvals)))

tickpref;

print_mfilename(mfilename)


