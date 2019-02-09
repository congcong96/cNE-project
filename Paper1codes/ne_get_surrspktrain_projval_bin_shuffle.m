function [spktrainPred, stamatPred, corrvals] = ne_get_surrspktrain_projval_bin_shuffle...
    (exp_site_nedata, sampleopt, nbins, stasig, sigopt, varargin)

% Shuffles spikes within projection value bins (15 or above) to get
% surrogate spiketrains predicted from a neuron's STA. The STA from the
% resulting surrogate spiketrain is validated by crosscorrelation: the STA
% from the 1st half of the actual spiketrain must be more correlated to the
% STA from the 2nd half of the surrogate spiketrain than that of the actual
% spiketrain (and vice versa, i.e. switching halves). 

% Updated 7/20/16 by JS. Automatic detection of sta_spiketrain for DF > 10

% Parse inputs and initialize some values
% ip = inputParser;
% addRequired(ip, 'exp_site_nedata', @isstruct)
% addOptional(ip, 'sampleopt', 'random', @ischar);
% addOptional(ip, 'nbins', 15, @isscalar);
% addOptional(ip, 'stasig', []);
% addOptional(ip, 'sigopt', 0, @isscalar);
% parse(ip, exp_site_nedata, sampleopt, nbins, stasig, sigopt)
% 
% exp_site_nedata = ip.Results.exp_site_nedata;
% sampleopt = ip.Results.sampleopt;
% nbins = ip.Results.nbins;
% stasig = ip.Results.stasig;
% sigopt = ip.Results.sigopt;


if nargin == 1
    sampleopt = 'random';
    nbins = 15;
    stasig = [];
    sigopt = 0;
elseif nargin == 2
    nbins = 15;
    stasig = [];
    sigopt = 0;
elseif nargin == 3
    stasig = [];
    sigopt = 0;
elseif nargin == 4
    sigopt = 0;
end

nedata = exp_site_nedata.nedata;
nf = nedata.nf;
nlags = nedata.nlags;
stim = regexp(exp_site_nedata.stim, 'rn\d{1,2}','match','once');

% initialize values for DFs > 10 and DFs <= 10
if exp_site_nedata.df > 10
    df = 10;
    spktrain_matrix = nedata.sta_spktrain;
else
    df = exp_site_nedata.df;
    spktrain_matrix = nedata.spktrain;
end
        
if isempty(varargin) && ~strcmp(sampleopt, 'coupratio')
    % get stimulus matrix
    drive = gcdr;
    subfolder = 'Ripple_Noise\downsampled_for_MID';
    stimfile = sprintf('%s-*10min_DFt1_DFf5_matrix.mat', stim);
    filename = dir(fullfile(drive,subfolder,stimfile));
    load(filename.name)

    % downsample
    stimulus = stim_mat(:,1:df:end);
    clear('stim_mat')

    % get stim matrix for STA computation
    stim = ne_create_stim_trial_from_stim_matrix(stimulus, spktrain_matrix, nlags);
elseif strcmp(sampleopt, 'coupratio')
    if length(varargin) == 1
        coupratio = varargin{1};
        drive = gcdr;
        subfolder = 'Ripple_Noise\downsampled';
        stimfile = sprintf('%s-*_DFt1_DFf5_matrix.mat', stim);
        filename = dir(fullfile(drive,subfolder,stimfile));
        load(filename.name)

        % downsample
        stimulus = stim_mat(:,1:df:end);
        clear('stim_mat')
        % get stim matrix for STA computation
        stim = ne_create_stim_trial_from_stim_matrix(stimulus, spktrain_matrix, nlags);
    else
        coupratio = varargin{1};
        stimulus = varargin{2};
        stim = varargin{3};
    end
else
    stimulus = varargin{1};
    stim = varargin{2};
end

% use normalized STAs if available
if isempty(stasig)
    stamat = nedata.stamat;
else
    stamat = stasig;
end

% get population firing rate if necessary
switch sampleopt
    case {'frdet','frprob','coupratio'}
        fr = mean(spktrain_matrix);
end

% initialize outputs
stamatPred = zeros(size(spktrain_matrix,1),size(stimulus,1)*nlags);
spktrainPred = zeros(size(spktrain_matrix));

for i = 1:size(spktrain_matrix, 1)
    fprintf('\nProcessing neuron %d of %d...', i, size(spktrain_matrix,1))
    
    locator = spktrain_matrix(i,:);
    sta = reshape(stamat(i,:), nf, nlags);

    % Calculate projection values
    [xprior, ~] = ca_calc_projection_from_sta_stimulus_spktrain(sta, stimulus, locator);
    
    % Get surrogate spiketrains 
    switch sampleopt
        case 'random'
            locatorPred = shufflePVbinsSpikePrediction(xprior, locator, nbins, 'random');
        case 'frdet'
            locatorPred = shufflePVbinsSpikePrediction(xprior, locator, nbins, 'frdet', fr);
        case 'frprob'
            locatorPred = shufflePVbinsSpikePrediction(xprior, locator, nbins, 'frprob', fr);
        case 'coupratio'
            locatorPred = couplingratioSpikePrediction (locator, xprior, fr, coupratio, nbins);
    end
    
    % Get surrogate spiketrains' STAs and their correlation with the actual
    % STAs
    staPred = locatorPred(nlags:end)*stim;
    
%     staPred = ca_calc_sta_from_stimulus_spktrain(stimulus, locatorPred, nlags);
%             
%     temp = corrcoef(sta(:), staPred(:));
%     corr(i) = temp(1,2);

    stamatPred(i,:) = staPred(:)';   
    spktrainPred(i,:) = locatorPred;

end

if sigopt == 1
    [stamatPred, ~] = ca_sig_sta_from_stim_obs_resp(stamatPred,...
        spktrainPred(:,nlags:end), stim, 100, 95);
end

corrvals = diag(corr(stamat', stamatPred'));

% downsample spiketrain to appropriate df for df > 10
if exp_site_nedata.df > 10
    dsfactor = exp_site_nedata.df/10;
    spktrainPred = downsample_spiketrain(spktrainPred, dsfactor);
    
    % remove last bin if incomplete
    if mod(size(spktrainPred,2), dsfactor) ~= 0
        spktrainPred = spktrainPred(:,1:end-1);
    end
    
end


fprintf('\n')
