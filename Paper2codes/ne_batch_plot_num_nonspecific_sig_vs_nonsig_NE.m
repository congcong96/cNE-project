function p = ne_batch_plot_num_nonspecific_sig_vs_nonsig_NE(nefiles)

% Haven't thought about what to do with this, but is potentially
% interesting.
% Non-specific events: cNE events that do not correspond with any of the
% cNE members spiking, i.e. driven mainly by non-specific population
% activity, probably during a huge UP state that doesn't involve cNE member
% neurons.

sigtab = cell(length(nefiles), 1);
nonsigtab = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    
    tab = ne_get_NE_event_member_spike_coincidence(exp_site_nedata);
    nonspecevents = cellfun(@(x) x{'0', {'Count'}}, tab);
    
    sigtab{i} = nonspecevents(exp_site_nedata.nedata.sig_NE_sta);
    nonsigtab{i} = nonspecevents(~exp_site_nedata.nedata.sig_NE_sta); 

end

nonspec.sig = cell2mat(sigtab);
nonspec.nonsig = cell2mat(nonsigtab);

p = plot_plotspread_and_boxplot(nonspec, 'ranksum/bonferroni');

