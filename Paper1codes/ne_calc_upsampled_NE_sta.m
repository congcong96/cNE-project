function NEsta = ne_calc_upsampled_NE_sta(actbin, stimtype)

stadf = 2;
nlags = 50;
% nf = 64;

drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);

stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stimtype);
stimfile = fullfile(stimfolder,stimfilename);
stimfile = dir(stimfile);
load(stimfile.name);

stimulus = stim_mat(:,1:stadf:end);
clear('stim_mat');

fn = fieldnames(actbin);

for i = 1:length(fn)
    
    fprintf('\nProcessing %d of %d bin sizes...', i, length(fn))
    
    if size(actbin.(fn{i}),2) == size(stimulus, 2)
        tempstim = stimulus;
    else %if size(sta_act,2) > size(sta_spktrain, 2)
        difference = size(stimulus,2) - size(actbin.(fn{i}),2);
        tempstim = stimulus(:, 1:end-difference);
    end
    
    stim = ne_create_stim_trial_from_stim_matrix(tempstim, [], nlags);    
    NEsta.(fn{i}) = actbin.(fn{i})(:,nlags:end) * stim;
    
    clear('stim')

end

fprintf('\n')
