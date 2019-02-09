function numpres = ne_plot_hist_of_unique_stims_per_repNE(repfiles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

numpres = cell(length(repfiles),1);

for i = 1:length(repfiles)
    
    load(repfiles{i})
    NEmem = repNE.NEmembers;
    NEidx = repNE.NEidx;
    
    numpres{i} = cellfun(@(x) length(unique(NEidx(x))), NEmem);    

end

presmat = categorical(cell2mat(numpres));
histogram(presmat, 'Normalization', 'probability');
xlim([0 3])
set(gca, 'xtick', 1:2, 'xticklabel', {'1','2'});
xlabel('Number of unique stimulus presentations represented in each NE')
ylabel('Ratio')
text(2, 0.7, sprintf('n = %d', length(presmat)));
tickpref;