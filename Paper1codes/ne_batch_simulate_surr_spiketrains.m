function surrspktrain = ne_batch_simulate_surr_spiketrains(exp_site_nedata, models, niter)

clc;
% change whenever necessary
coupratio = 1;

nedata = exp_site_nedata.nedata;
df = 10; % temp for ne calculation only
nlags = nedata.nlags;
stimtype = regexp(exp_site_nedata.stim, 'rn\d{1,2}','match','once');

drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfile = sprintf('%s-*_DFt1_DFf5_matrix.mat', stimtype);
filename = dir(fullfile(drive,subfolder,stimfile));
load(filename.name)

% downsample
stimulus = stim_mat(:,1:df:end);
clear('stim_mat')

% get stim matrix for STA computation
stim = ne_create_stim_trial_from_stim_matrix(stimulus, [], nlags);

for i = 1:length(models)
    fprintf('Simulating %s...\n', models{i})
    
    model = models{i};
    spktrainPred = cell(niter, 1);
    for j = 1:niter
        fprintf('Iteration %d of %d\n', j, niter)
    
        switch model
            case 'preISI' % preserved ISI
                spktrainPred{j} = ne_circularly_shuffle_spkmatrix (exp_site_nedata);
            case 'preRF' % preserved RF
                [spktrainPred{j}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'random', 15, [], 0, stimulus, stim);
            case 'preRFnorm' % preserved RF with RF normalization
                [spktrainPred{j}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'random', 15, stamat);
            case 'prePFRRF' % preserved PFR + RF
                [spktrainPred{j}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'frprob', 15, [], 0, stimulus, stim);
            case 'prePFRRFnorm' % preserved PFR + RF with RF normalization
                [spktrainPred{j}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'frprob', 15, stamat);        
            case 'prePFR' % preserved PFR
                spktrainPred{j} = ca_get_surrspktrain_weighted_popfr_shuffle (exp_site_nedata, 1);
            case 'jitter' % jitter method
                spktrainPred{j} = ca_get_jittered_spktrains(exp_site_nedata);
            case 'coupratio'
                spktrainPred{j} = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'coupratio', 15, [], 0, coupratio, stimulus, stim);
            otherwise
                error('Inappropriate model option')
        end
        
    end
    
    surrspktrain.(sprintf('%s',model)) = spktrainPred;
    clc;
end