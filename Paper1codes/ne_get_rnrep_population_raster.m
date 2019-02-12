function corrvals = get_rnrep_population_raster(spk, trigger, binsize, plotopt, neuronopt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 4
    plotopt = 0;
end

offset = 200;
triallength = 4700;
numtrials = 50;

if exist('neuronopt', 'var') && ~isempty(neuronopt)
    spk = spk(neuronopt);
end

binedges = 0:binsize:triallength;

trialspkbin = zeros(numtrials, length(binedges) - 1, length(spk));

for i = 1:length(spk)
    
    trialspk = ne_calc_rnrep_trialspk(spk(i).spiketimes, trigger, offset);
    trialspkbin(:,:,i) = cell2mat(cellfun(@(x) histcounts(x, binedges),...
        trialspk, 'UniformOutput', 0)');
    
end

popbin = sum(trialspkbin,3);
    
if plotopt == 1
    figure;
    imagesc(popbin)
    colormap 'jet'
end
% switch neopt
%     case 'raw'
%         popbin = reshape(zscore(popbin(:)), size(popbin,1), size(popbin,2));
%         if plotopt == 1
%             imagesc(popbin)
%             colormap 'jet'
%         end
%     case 'logical'
%         [~, idx] = sort(popbin(:), 'descend');
%         popbin = zeros(size(popbin));
%         popbin(idx(1:numspks)) = 1;
%         plotSpikeRaster(logical(popbin),'PlotType','vertline2', 'TimePerBin', binsize/1000);
%         
% end

temp = triu(corr(popbin'),1);
corrvals = temp(temp ~=0);   

