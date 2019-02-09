function [shuffledspkmatrix, varargout] = ne_shuffle_rnrep_trial_order(spk, trigger, binsize)

offset = 200;
triallength = 4700;

shuffledspkmatrix = zeros(length(spk), triallength/binsize*length(trigger));
realspkmatrix = zeros(length(spk), triallength/binsize*length(trigger));
randidx = zeros(length(spk),length(trigger));

for i = 1:length(spk)
    
    trialspk = ne_calc_rnrep_trialspk(spk(i).spiketimes, trigger, offset);
    trialspkbin = cell2mat(cellfun(@(x) histcounts(x, 0:binsize:triallength),...
        trialspk, 'UniformOutput', 0)');
    randidx(i,:) = randperm(length(trigger));
    shufftrialspkbin = trialspkbin(randidx(i,:),:)';
    shuffledspkmatrix(i,:) = shufftrialspkbin(:)';
    trialspkbin = trialspkbin';
    realspkmatrix(i,:) = trialspkbin(:)';
    
end

varargout{1} = randidx;
varargout{2} = realspkmatrix;

return