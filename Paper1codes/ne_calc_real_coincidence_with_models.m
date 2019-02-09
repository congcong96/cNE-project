function [neurons, ratio, varargout] = ne_calc_real_coincidence_with_models...
    (exp_site_nedata, combopt, numsample, modelopt, niter, NEopt, coupratio, varargin)
% Updated 7/21/16 by JS

%% Parse input arguments and initialize nedata

% p = inputParser;
% addOptional(p, 'combopt', 'sameNE', @isstr);
% addOptional(p, 'numsample', 100, @isscalar);
% addOptional(p, 'modelopt', 'coupratio', @isstr);
% addOptional(p, 'niter', 1, @isscalar);
% addOptional(p, 'NEopt', 1, @isscalar);
% addOptional(p, 'coupratio', 0.5, @(x)validateattributes(x, {'double'},{'>=', 0, '<=', 1}));
% parse(p, combopt, numsample, modelopt, niter, NEopt, coupratio);

if nargin == 2
    numsample = 100;
    modelopt = 'coupratio';
    niter = 1;
    NEopt = 1;
    coupratio = 0.5;
elseif nargin == 3
    modelopt = 'coupratio';
    niter = 1;
    NEopt = 1;
    coupratio = 0.5;
elseif nargin == 4
    niter = 1;
    NEopt = 1;
    if strcmp(modelopt, 'coupratio')
        coupratio = 0.5;
    else
        coupratio = [];
    end
elseif nargin == 5
    NEopt = 1;
    if strcmp(modelopt, 'coupratio')
        coupratio = 0.5;
    else
        coupratio = [];
    end
elseif nargin == 6
    if strcmp(modelopt, 'coupratio')
        coupratio = 0.5;
    else
        coupratio = [];
    end
end


if ischar(modelopt)
    if strcmp(modelopt, 'all')
        modelopt = {'repeat','preISI','preRF','prePFR','prePFRRF','jitter'};
    elseif strcmp(modelopt, 'paper')
        modelopt = {'repeat','preISI','preRF','prePFRRF'};
    else
        modelopt = {modelopt};
    end
end


nedata = exp_site_nedata.nedata;


%% Get combinations of neurons and initialize outputs

% if combinations of neurons are specified
if iscell(combopt)   
    
    neurons = combopt;
    
else
    
% else get combinations of neurons in the same NE or otherwise

    numcomp = 2:5; %number of neurons to consider at once
    neurons = cell(length(numcomp),1);
    
    switch combopt
        
        case 'sameNE'
            % get NE combinations
            for i = 1:length(numcomp)
                neurons{i} = ne_find_NE_pairs_or_groups(exp_site_nedata, numcomp(i), numsample, 1);       
            end
            
    
        case 'diffNE'
            % get nonNE combinations
            for i = 1:length(numcomp)
                [neurons{i},~] = ne_find_non_NE_pairs_or_groups(exp_site_nedata, numcomp(i), numsample, 1);
            end
        otherwise
            error('Invalid combination option.')
    end
    
end

%% Get coincidence ratios / NE stats (if needed) for real spiketrain
rspktrain = nedata.spktrain;
ratio.real_ratio = cell(length(neurons),1);

for ii = 1:length(neurons)
    ratio.real_ratio{ii} = ne_calc_coincidence_within_spktrain(rspktrain, neurons{ii});
end

if NEopt == 1
    nedata = ca_spkmatrix_to_ensembles(rspktrain);
    NEstats.real.numNE = size(nedata.ensembles,2);
    NEstats.real.numsig = sum(nedata.eigenvalues...
        > nedata.lambda_max | nedata.eigenvalues < ...
        nedata.lambda_min);
end

%% Normalize STAs if needed
if sum(ismember(modelopt,{'preRFnorm','prePFRRFnorm'})) > 0

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
    [stamat, ~] = ca_sig_sta_from_stim_obs_resp(stamat, resp, stim, 100, 95);
end
    

%% Get surrogate spiketrains
% for all models
for j = 1:length(modelopt)
    model = modelopt{j};
        
    switch model
        case 'repeat' % for repeated stimuli which does not require iterations
            surrratio = cell(length(neurons), 1);
            spktrainPred = cell(length(varargin),1);
            for iii = 1:length(varargin)                    
                spktrainPred{iii} = varargin{iii}.nedata.spktrain;
            end
                
            for k = 1:length(neurons)
                
                surrratio{k} = ne_calc_coincidence_across_spktrains(rspktrain, spktrainPred, neurons{k});
                
            end
                
        otherwise % for all other models requiring iterations
            % for all iterations
            surrratio = cell(length(neurons), niter);
            spktrainPred = cell(niter,1);
            
            for jjj = 1:niter  
                % get surrogate spiketrains with different models
                switch model
                    case 'preISI' % preserved ISI
                        spktrainPred{jjj} = ne_circularly_shuffle_spkmatrix (exp_site_nedata);
                    case 'preRF' % preserved RF
                        [spktrainPred{jjj}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'random');
                    case 'preRFnorm' % preserved RF with RF normalization
                        [spktrainPred{jjj}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'random', 15, stamat);
                    case 'prePFRRF' % preserved PFR + RF
                        [spktrainPred{jjj}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'frprob');
                    case 'prePFRRFnorm' % preserved PFR + RF with RF normalization
                        [spktrainPred{jjj}, ~,~] = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'frprob', 15, stamat);        
                    case 'prePFR' % preserved PFR
                        spktrainPred{jjj} = ca_get_surrspktrain_weighted_popfr_shuffle (exp_site_nedata, 1);
                    case 'jitter' % jitter method
                        spktrainPred{jjj} = ca_get_jittered_spktrains(exp_site_nedata);
                    case 'coupratio'
                        spktrainPred{jjj} = ne_get_surrspktrain_projval_bin_shuffle (exp_site_nedata, 'coupratio', 15, [], 0, coupratio);
                    otherwise
                        error('Inappropriate model option')
                end
                
                %% Calculate NE stats if required
                if NEopt == 1
                    nedata = ca_spkmatrix_to_ensembles(spktrainPred{jjj});
                    NEstats.(sprintf('%s',model)).numNE(jjj) = size(nedata.ensembles,2);
                    NEstats.(sprintf('%s',model)).numsig(jjj) = ...
                        sum(nedata.eigenvalues > nedata.lambda_max |...
                        nedata.eigenvalues < nedata.lambda_min);
                end
                        


%% Coincidence ratio calculations for modelled spiketrains

                for k = 1:length(neurons)                    
                    surrratio{k, jjj} = ne_calc_coincidence_within_spktrain(spktrainPred{jjj}, neurons{k});
                end


            end
    end
    
    if NEopt == 1
        varargout{1} = NEstats;
        varargout{2}.(sprintf('%s_spktrain',model)) = spktrainPred;
    else
        varargout{1}.(sprintf('%s_spktrain',model)) = spktrainPred;
    end
    
    ratio.(sprintf('%s_ratio',model)) = surrratio;
    
end

    
end  


