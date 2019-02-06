function ne_plot_pca_waveforms_scatter(waveforms, idx, label, color)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(idx)
    proj = waveforms.projections(idx);
    cat = waveforms.category(idx);
else
    proj = waveforms.projections;
    cat = waveforms.category;
end

maxuniq = cellfun(@(x) max(length(x)), cat);

if isempty(color)
    color = distinguishable_colors(maxuniq);
end

figure;
c = 0;

for i = 1:length(proj)
    
    c = c + 1;
    if c == 5
        c = 1;
        figure;
    end
    
    subplot(2,2,c)
    hold on
    
    uniqcat = unique(cat{i});
    
    if isempty(color)
        color = distinguishable_colors(length(uniqcat));
    end
        
    
    for j = 1:length(uniqcat)
        scatter(proj{i}(cat{i}==uniqcat(j), 1), proj{i}(cat{i}==uniqcat(j), 2), 6, color(j,:), 'filled');
    end
       
    if ~isempty(label)
        legend(label(1:length(uniqcat)))
    end 
    

end

