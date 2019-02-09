function bicell = ne_get_bicellular_spktrain_subset(exp_site_nedata, spktrain)



%downsample spktrain to 5dft
dsspktrain = downsample_spiketrain(spktrain, 10);

bispktrain = get_bicellular_spktrains(dsspktrain);
bicomb = cell2mat({bispktrain.comb}');

NEneurons = ne_find_NE_pairs_or_groups(exp_site_nedata, 2);

[~,NEidx,~] = intersect(bicomb,NEneurons,'rows');
nonNEidx = setdiff(1:size(bicomb,1), NEidx);

bicell.NE = bispktrain(NEidx);
bicell.nonNE = bispktrain(nonNEidx);

% for i = 1:size(dsspktrain,1)
%     
%     NEidx = logical(cellfun(@(x) sum(i == x), NEmembers));
%     NEneurons = setdiff(cell2mat(NEmembers(NEidx)), i);
%     nonNEneurons = setdiff(1:size(dsspktrain,1), [NEneurons; i]);
%     
%     
%     
% end
    