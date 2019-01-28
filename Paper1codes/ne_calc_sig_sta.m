function sig_stamat = ne_calc_sig_sta(exp_site_nedata)

nedata = exp_site_nedata.nedata;
stamat = nedata.stamat;
if exp_site_nedata.df > 10
    df = 10;
else
    df = exp_site_nedata.df;
end

stimtype = regexp(exp_site_nedata.stim, 'rn\d{1,2}','match','once');
drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfile = sprintf('%s-*_DFt1_DFf5_matrix.mat', stimtype);
filename = dir(fullfile(drive,subfolder,stimfile));
stimmat = load(filename.name);
stim_mat = stimmat.stim_mat;
stimulus = stim_mat(:,1:df:end);
clear('stim_mat')
if isfield(nedata, 'sta_spktrain')
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(stimulus, nedata.sta_spktrain, nedata.nlags);
else
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(stimulus, rspktrain, nedata.nlags);
end
[sig_stamat, ~] = ca_sig_sta_from_stim_obs_resp(stamat, resp, stim, 100, 95);
