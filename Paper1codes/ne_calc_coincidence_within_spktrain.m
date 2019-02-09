function coinratio = ne_calc_coincidence_within_spktrain(spktrain, comb)

% Calculates coincidence ratio based on spiketrain and a matrix of
% combinations, where each row represents one combination

% Written 7/21/16 by JS

% initialize output
coinratio = zeros(size(comb,1),1);
% get number of neurons per combination
num = size(comb,2);

for i = 1:size(comb,1)
    
    % get spiketrains of neurons in one combination
    compspktrain = spktrain(comb(i,:),:);
    % get numerator and denominator of coincidence ratio
    nume = sum(sum(compspktrain >= 1) == num);
    deno = min(sum(compspktrain >= 1, 2));
    % coincidence ratio
    coinratio(i) = nume/deno;

end

