function actstruct = ne_interpolate_upsampled_activities_for_diff_binsizes(NEact, stadf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('stadf', 'var')
    stadf = 2;
end
% nlags = 50;
% nf = 64;

% drive = gcdr;
% subfolder = 'Ripple_Noise\downsampled';
% stimfolder = fullfile(drive,subfolder);
% 
% stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stimtype);
% stimfile = fullfile(stimfolder,stimfilename);
% stimfile = dir(stimfile);
% load(stimfile.name);
% 
% stimulus = stim_mat(:,1:stadf:end);
% clear('stim_mat');

% [stim] = ne_create_stim_trial_from_stim_matrix(stimulus, [], nlags);

fn = fieldnames(NEact);

for i = 1:length(fn)
    
    fprintf('\nProcessing %d of %d bin sizes...', i, length(fn))
    
    activity = NEact.(fn{i});
    df = str2double(regexp(fn{i}, '(?<=(^df))\d{1,3}','match','once'));
    
    for j = 1:size(activity,1)
        actstruct.(fn{i})(j,:) = interp(activity(j,:), df/stadf);
    end
    
%     if size(act_cell{i},2) == size(stimulus, 2)
%         NEsta.(fn{i}) = act_cell{i}(:,nlags:end)* stim;
%     else %if size(sta_act,2) > size(sta_spktrain, 2)
%         difference = size(stimulus,2) - size(act_cell{i},2);
%         NEsta.(fn{i}) = act_cell{i}(:,nlags:end)* stim(1:end-difference,:);
%     end

end



% act_mat = cell2mat(act_cell);
% 
% temp = act_mat(:,nlags:end) * stim;
% clear('act_mat','stim');
% 
% idxend = cumsum(cellfun(@(x) size(x,1), act_cell));
% idxstart = [1 ;idxend(1:end-1)+1];
% 
% for k = 1:length(fn)
%     NEsta.(fn{k}) = temp(idxstart(k):idxend(k),:);
% end
fprintf('\n')


