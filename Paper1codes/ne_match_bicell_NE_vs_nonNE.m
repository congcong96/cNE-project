function matchedbicell = ne_match_bicell_NE_vs_nonNE(bicell)

matchedbicell(length(bicell.NE)).NEspktrain = [];

% nonNEcomb = {bicell.nonNE.comb}';
NEcomb = {bicell.NE.comb}';

c = 1;

for i = 1:length(bicell.nonNE)
    match1 = cellfun(@(x) any(bicell.nonNE(i).comb(1) == x), NEcomb);
    match2 = cellfun(@(x) any(bicell.nonNE(i).comb(2) == x), NEcomb);
        
    idx1 = find(match1 == 1);
    idx2 = find(match2 == 1);
    idx = [idx1; idx2];
    
    if isempty(idx)
        
        continue
        
    elseif length(idx) == 1
        chosen = idx;
    else
        chosen = randsample(idx,1);
    end
    
    
%     if sum(match1) >= sum(match2)
%         idx = find(match1 == 1);
%         if length(idx) == 1
%             chosen = idx;
%         else
%             chosen = randsample(idx, 1);
%         end
%     else
%         idx = find(match2 == 1);
%         if length(idx) == 1
%             chosen = idx;
%         else
%             chosen = randsample(idx, 1);
%         end
%     end
    
    matchedbicell(c).nonNEspktrain = bicell.nonNE(i).spktrain;
    matchedbicell(c).nonNEcomb = bicell.nonNE(i).comb;
    matchedbicell(c).NEspktrain = bicell.NE(chosen).spktrain;
    matchedbicell(c).NEcomb = bicell.NE(chosen).comb;
    matchedbicell(c).refneuron = intersect(matchedbicell(c).NEcomb, matchedbicell(c).nonNEcomb);
    matchedbicell(c).nonNEneuron = setdiff(matchedbicell(c).nonNEcomb, matchedbicell(c).refneuron);
    matchedbicell(c).NEneuron = setdiff(matchedbicell(c).NEcomb, matchedbicell(c).refneuron);
    
    %update comb
%     nonNEcomb(chosen) = [];
%     bicell.nonNE(chosen) = [];
    c = c+1;
end

    