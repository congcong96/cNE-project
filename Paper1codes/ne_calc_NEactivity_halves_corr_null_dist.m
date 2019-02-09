function [actcorr, shuffledactcorr] = ne_calc_NEactivity_halves_corr_null_dist(exp_site_nedata)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


[actcorr] = ne_predict_subset_of_NE_activity_with_ICweights(exp_site_nedata);

nedata = exp_site_nedata.nedata;
spktrain = nedata.spktrain;

% requiredsamp = 10000;
% 
% %compute min number of iterations required to get at least 10000 samples
% corrvalspermat = size(actcorr,1)^2;
% niter = ceil(requiredsamp / corrvalspermat);
% 
% shuffledactcell = ne_calc_shuffled_spktrain_halves_NEact(spktrain,'circular',niter,'splitopt', 'contiguous', 'actopt', 1);
% close all
% 
% shuffledactcorr = cell2mat(shuffledactcell);
% shuffledactcorr = shuffledactcorr(:);
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

