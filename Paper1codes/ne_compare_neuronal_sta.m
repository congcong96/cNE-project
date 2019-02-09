function [autocorrs, crosscorrs] = ne_compare_neuronal_sta(spktrain, stim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

spktrainhalf1 = spktrain(:,1:floor(size(spktrain,2)/2));
spktrainhalf2 = spktrain(:,floor(size(spktrain,2)/2)+1:end);

nlags = 50;
stadf = 2;

% get stim_mat for calculating STA (long stim)
currdrive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimmatfile = gfn([fullfile(currdrive, subfolder) sprintf('\\%s-*DFt1_DFf5_matrix.mat', stim)]);
load(stimmatfile{1})
stim_mat = stim_mat(:,1:stadf:end);
stimhalf1 = stim_mat(:,1:floor(size(stim_mat,2)/2));
stimhalf2 = stim_mat(:,floor(size(stim_mat,2)/2)+1:end);

% if size(NEactproj,2) ~= size(NEactUS,2)
%     difference = size(NEactproj,2) - size(NEactUS,2);
%     NEactproj = NEactproj(:,1:end-difference);
%     stim_mat = stim_mat(:,1:end-difference);
% end

stahalf1 = quick_calc_sta(stimhalf1, spktrainhalf1, nlags);
stahalf2 = quick_calc_sta(stimhalf2, spktrainhalf2, nlags);

% [sigstahalf1, ~] = ca_sig_sta_from_stim_obs_resp(stahalf1, NEacthalf1, stimhalf1, 20, 95);
% [sigstahalf2, ~] = ca_sig_sta_from_stim_obs_resp(stahalf2, NEacthalf2, stimhalf2, 20, 95);



% staUS = quick_calc_sta(stim_mat, NEactUS, nlags);

% figure; hold on; plot(NEactUS(1,:)); plot(NEactproj(1,:));

% normstacorr = corr(sigstahalf1', sigstahalf2');
stacorr = corr(stahalf1', stahalf2');

for i = 1:length(stacorr)
    r = stacorr(i,:);
    c = stacorr(:,i);
    r(i) = [];
    c(i) = [];
    crosscorrs{i} = [c;r'];
    realval = repmat(stacorr(i,i),(length(stacorr)-1)*2, 1);
    [~,p(i,1)] = ttest(realval, crosscorrs{i});
%     p(i,1) = signrank(realval,a);
    
%     r = normstacorr(i,:);
%     c = normstacorr(:,i);
%     r(i) = [];
%     c(i) = [];
%     a = [c;r'];
%     realval = repmat(normstacorr(i,i),(length(normstacorr)-1)*2, 1);
% %     [~,p(i,2)] = ttest(realval, a);
%     p(i,2) = signrank(realval,a);
        

end

autocorrs = diag(stacorr);


end
