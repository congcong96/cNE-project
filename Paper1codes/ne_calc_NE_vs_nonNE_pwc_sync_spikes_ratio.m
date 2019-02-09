function [NEsyncratio, nonNEsyncratio] = ne_calc_NE_vs_nonNE_pwc_sync_spikes_ratio(exp_site_nedata, allpwc)


% get all NE pairs
NEneurons = ne_find_NE_pairs_or_groups(exp_site_nedata, 2);

% calc sync spikes
syncspikes = calc_pwc_synchronous_spikes_ratio(allpwc, 10);
pairs = cell2mat({syncspikes.pairs}');

[~,~,NEidx] = intersect(NEneurons,pairs,'rows');
nonNEidx = setdiff(1:length(syncspikes), NEidx);
NEsyncratio = [syncspikes(NEidx).ratio];
nonNEsyncratio = [syncspikes(nonNEidx).ratio];

end

