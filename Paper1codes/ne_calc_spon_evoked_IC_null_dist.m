function [shuffledcorrvals, sponcorrvals, evokedcorrvals] = ne_calc_spon_evoked_IC_null_dist(spon_nedata, evoked_nedata)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

spon_spktrain = spon_nedata.spktrain;
evoked_spktrain = evoked_nedata.spktrain;

requiredsamp = 100000;

%compute min number of iterations required to get at least 100000 samples
corrvalspermat = size(spon_nedata.Patterns,2) * size(evoked_nedata.Patterns,2);
numcomb = requiredsamp/corrvalspermat;
niter = ceil(sqrt(numcomb));



sponICshuffle = ne_calc_shuffled_spktrain_ICs(spon_spktrain,'circular', niter);
evokedICshuffle = ne_calc_shuffled_spktrain_ICs(evoked_spktrain, 'circular', niter);
close all

shuffledcorrvals = [];

for i = 1:niter
    
    for j = 1:niter
        
        temp = corr(squeeze(sponICshuffle(:,:,i)), squeeze(evokedICshuffle(:,:,j)));
        shuffledcorrvals = [shuffledcorrvals; temp(:)];
        
    end
end

sponevokedcorrmat = corr(spon_nedata.Patterns, evoked_nedata.Patterns);

    
[~,idx] = max(abs(sponevokedcorrmat),[],2);
sponcorrvals = cellfun(@(x,y) sponevokedcorrmat(x,y), num2cell(1:size(sponevokedcorrmat,1)), num2cell(idx'));

[~,idx] = max(abs(sponevokedcorrmat),[],1);
evokedcorrvals = cellfun(@(x,y) sponevokedcorrmat(x,y), num2cell(idx), num2cell(1:size(sponevokedcorrmat,2)));



end

