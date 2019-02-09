function [coinratio, varargout] = ne_calc_coincidence_across_spktrains(maintrain, othertrains, comb, compareopt, comblim, mtype)

% Calculates coincidence ratio based on a few spiketrains and a matrix of
% combinations, where each row represents one combination, and where each
% neuron in each combination comes from a different spiketrain.

% Written 7/21/16 by JS, updated 8/30/16

if nargin == 3
    compareopt = 'fixedcomb';
    comblim = [];
    mtype = [];
end

if size(comb,2) > length(othertrains) + 1
    error('Too few spiketrains for the number of neurons per combination.')
end

% initialize output
coinratio = zeros(size(comb,1),1);

% get number of combinations and number of neurons per combination
numcompare = size(comb,2);
numcomb = size(comb,1);

switch compareopt
    case 'allcomb'
        numcombi = nchoosek(1:length(othertrains), numcompare-1);
        if size(numcombi,1) < comblim
            temp = [zeros(size(numcombi,1),1) numcombi];
            traincomb = [];
            for i = 1:size(temp,1)
                traincomb = [traincomb; perms(temp(i,:))];
            end
            if size(traincomb,1) > comblim
                randidx = randsample(size(traincomb,1), comblim);
                traincomb = traincomb(randidx,:);
            end
        else %size(numcombi,1) >= comblim
           randidx = randsample(size(numcombi,1), comblim);
           traincomb = [zeros(comblim,1) numcombi(randidx,:)];
           for i = 1:comblim
               permidx = randperm(numcompare);
               traincomb(i,:) = traincomb(i,permidx);
           end
        end
end

for k = 1:numcomb
    switch compareopt
        case 'fixedcomb'
           
            % keeps identity of neurons, i.e. 2nd neuron from 2nd spktrain,
            % etc.
            mainneuron = randsample(comb(k,:),1);
            compspktrain(1,:) = maintrain(mainneuron,:);
            
            otherneurons = comb(k,:);
            otherneurons(otherneurons == mainneuron) = [];
            otherneurons = otherneurons(randperm(numcompare-1));
            
            for j = 1:length(otherneurons)
                compspktrain(j+1,:) = othertrains{j}(otherneurons(j),:);
            end
            
            % get numerator and denominator of coincidence ratio
            nume = sum(sum(compspktrain >= 1) == numcompare);
            deno = min(sum(compspktrain >= 1, 2));
            % coincidence ratio
            coinratio(k) = nume/deno;
            
            if isnan(coinratio(k))
                coinratio(k) = 0;
            end
            
        case 'onecomb'
            % randomly choose one main neuron from main spiketrain
            mainneuron = randsample(comb(k,:),1);
            compspktrain(1,:) = maintrain(mainneuron,:);
            
            % assign other neurons from other spiketrains randomly
            otherneurons = comb(k,:);
            otherneurons(otherneurons == mainneuron) = [];
            otrainorder = randsample(length(othertrains), numcompare - 1);
                        
            for j = 1:length(otherneurons)
                compspktrain(j+1,:) = othertrains{otrainorder(j)}(otherneurons(j),:);
            end
            
            % get numerator and denominator of coincidence ratio
            nume = sum(sum(compspktrain >= 1) == numcompare);
            deno = min(sum(compspktrain >= 1, 2));
            % coincidence ratio
            coinratio(k) = nume/deno;

        case 'allcomb' % get all permutations of neurons (up to comb limit) per combination
            fprintf('\nProcessing combination %d of %d...',k,numcomb)
           for j = 1:size(traincomb,1)
               for ii = 1:size(traincomb,2)
                   if traincomb(j,ii) == 0
                       compspktrain(ii,:) = maintrain(comb(k,ii),:);
                   else
                       compspktrain(ii,:) = othertrains{traincomb(j,ii)}(comb(k, ii),:);
                   end
               end
                % get numerator and denominator of coincidence ratio
                nume = sum(sum(compspktrain >= 1) == numcompare);
                deno = min(sum(compspktrain >= 1, 2));
                allcoinratio(j,k) = nume/deno;
           end
           switch mtype
               case 'median'
                   coinratio(k) = median(allcoinratio(:,k));
               case 'mean'
                   coinratio(k) = mean(allcoinratio(:,k));
           end
            
    end
    
%     % get numerator and denominator of coincidence ratio
%     nume = sum(sum(compspktrain >= 1) == numcompare);
%     deno = min(sum(compspktrain >= 1, 2));
%     % coincidence ratio
%     coinratio(i) = nume/deno;
    
end

if strcmp(compareopt, 'allcomb')
    varargout{1} = allcoinratio;
    fprintf('\n')
end
            