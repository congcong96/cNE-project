function [corrmat, shuffledcorrvals] = ne_calc_IC_halves_null_dist(exp_site_nedata)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


[corrmat, IC1, IC2] = ne_compare_ICweights_of_halves(exp_site_nedata, 'interleaved');

nedata = exp_site_nedata.nedata;
spktrain = nedata.spktrain;

requiredsamp = 10000;

%compute min number of iterations required to get at least 10000 samples
corrvalspermat = size(IC1,2) * size(IC2,2);
niter = ceil(requiredsamp / corrvalspermat);

ICshufflecell = ne_calc_shuffled_spktrain_halves_ICs(spktrain,'circular',niter);
close all

shuffledcorrvals = [];

for i = 1:niter

    temp = corr(ICshufflecell{i,1},ICshufflecell{i,2});
    shuffledcorrvals = [shuffledcorrvals; temp(:)];
            
end

% sponevokedcorrmat = corr(nedata.spon_patterns, nedata.evoked_patterns);
% 
%     
% [~,idx] = max(abs(sponevokedcorrmat),[],2);
% sponcorrvals = cellfun(@(x,y) sponevokedcorrmat(x,y), num2cell(1:size(sponevokedcorrmat,1)), num2cell(idx'));
% 
% [~,idx] = max(abs(sponevokedcorrmat),[],1);
% evokedcorrvals = cellfun(@(x,y) sponevokedcorrmat(x,y), num2cell(idx), num2cell(1:size(sponevokedcorrmat,2)));
% 


end

