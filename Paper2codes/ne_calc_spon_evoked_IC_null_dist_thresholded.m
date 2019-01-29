function [shuffledcorrvals, sponcorrvals, evokedcorrvals] = ne_calc_spon_evoked_IC_null_dist_thresholded(spon_nedata, evoked_nedata)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

CI = evoked_nedata.CI;
spon_ICs = spon_nedata.Patterns <= CI(1) | spon_nedata.Patterns >= CI(2);
evoked_ICs = evoked_nedata.Patterns <= CI(1) | evoked_nedata.Patterns >= CI(2);

requiredsamp = 100000;

%compute min number of iterations required to get at least 100000 samples
corrvalspermat = size(spon_ICs,2) * size(evoked_ICs,2);
numcomb = requiredsamp/corrvalspermat;
niter = ceil(sqrt(numcomb));

sponICshuffle = zeros([size(spon_ICs) niter]);
evokedICshuffle = zeros([size(evoked_ICs) niter]);

for i = 1:niter
    
    for j = 1:size(spon_ICs, 2)
        idx = randsample(size(spon_ICs, 1), size(spon_ICs, 1));
        sponICshuffle(:,j,i) = spon_ICs(idx, j);
    end
    
    for j = 1:size(evoked_ICs, 2)
        idx = randsample(size(evoked_ICs, 1), size(evoked_ICs, 1));
        evokedICshuffle(:,j,i) = evoked_ICs(idx, j);
    end
end

shuffledcorrvals = zeros(corrvalspermat * niter^2, 1);
c = 1;

for i = 1:niter
    
    for j = 1:niter
        
        temp = corr(squeeze(sponICshuffle(:,:,i)), squeeze(evokedICshuffle(:,:,j)));
        shuffledcorrvals((c-1)*corrvalspermat+1:c*corrvalspermat) = temp(:);
        c = c+1;
        
    end
end

sponevokedcorrmat = corr(spon_nedata.Patterns, evoked_nedata.Patterns);

    
[~,idx] = max(abs(sponevokedcorrmat),[],2);
sponcorrvals = cellfun(@(x,y) sponevokedcorrmat(x,y), num2cell(1:size(sponevokedcorrmat,1)), num2cell(idx'));

[~,idx] = max(abs(sponevokedcorrmat),[],1);
evokedcorrvals = cellfun(@(x,y) sponevokedcorrmat(x,y), num2cell(idx), num2cell(1:size(sponevokedcorrmat,2)));



end

