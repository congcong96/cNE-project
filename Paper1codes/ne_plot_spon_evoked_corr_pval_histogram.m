function ne_plot_spon_evoked_corr_pval_histogram(corrvalstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

edges = 0:0.05:1;

sponp = cell2mat({corrvalstruct.sponpvals});
evokedp = cell2mat({corrvalstruct.evokedpvals});


figure;
histogram(sponp, edges, 'Normalization', 'probability');
xlabel('P-value')
ylabel('Ratio')
text(0.5, 0.5, sprintf('n = %d', length(sponp)))
title('Spontaneous correlations p-values')
tickpref;

figure;
histogram(evokedp, edges, 'Normalization', 'probability');
xlabel('P-value')
ylabel('Ratio')
text(0.5, 0.5, sprintf('n = %d', length(evokedp)))
title('Evoked correlations p-values')
tickpref;


end

