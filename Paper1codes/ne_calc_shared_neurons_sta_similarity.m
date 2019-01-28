function corrstruct = ne_calc_shared_neurons_sta_similarity (exp_site_nedata, varargin)

% Gets null distribution of Pearson's correlation values between random
% subsets of spikes for each neuron that is present in more than one cNE,
% and actual Pearson's correlation values between subsets of spikes based
% on cNEs. Number of spikes sampled depends on hte lowest number of spikes
% for each neuron in each cNE.

% staopt: If 'raw', do not threshold STA before comparing. If 'norm',
%   threshold STA before comparing.
% compopt: If 'pairs', compare only non-overlapping sets of spikes. If
%   'anycomb', make 100 random subset spike trains and do all pairwise 
%   correlations.
% sigopt: If true, calculate only for neurons with significant STAs.
% threshalpha: Determine threshold for 'binarizing' cNE train.
% spkthresh: Only calculate STA similarity if spike count of smallest
%   subset is larger than this threshold.

% Updated 5/19/17 by JS
% Updated 7/23/18 by JS, updated to use with quick_calc_sta.m, updated to
% add spkthresh option in ne_get_shared_neuron_spktrain_subsets.m and
% updated description.

ip = inputParser;
addRequired(ip, 'exp_site_nedata', @isstruct)
addParameter(ip, 'staopt', 'raw',  @(x) strcmp(x, 'raw') || strcmp(x, 'norm'))
addParameter(ip, 'compopt', 'pairs', @(x) strcmp(x, 'pairs') || strcmp(x, 'anycomb'))
addParameter(ip, 'sigopt', 1, @(x) x == 0 || x == 1)
addParameter(ip, 'threshalpha', 99.5, @(x) x >= 99 && x <= 99.9)
addParameter(ip, 'spkthresh', 100, @isscalar);
parse(ip, exp_site_nedata, varargin{:})

exp_site_nedata = ip.Results.exp_site_nedata;
staopt = ip.Results.staopt;
compopt = ip.Results.compopt;
sigopt = ip.Results.sigopt;
threshalpha = ip.Results.threshalpha;
spkthresh = ip.Results.spkthresh;

spkstruct = ne_get_shared_neuron_spktrain_subsets(exp_site_nedata, 'threshalpha', threshalpha, 'spkthresh', spkthresh);

if isempty(spkstruct)
    corrstruct = [];
    return
end

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;

if sigopt
    sig_neurons = nedata.sig_neuron_sta([spkstruct.neuron]);
    spkstruct = spkstruct(sig_neurons);
    
    if isempty(spkstruct)
        warning('Dataset has no shared neurons with significant STAs.')
        corrstruct = [];
        return
    end
end
    

switch compopt %option to only compare pairs of non-overlapping spikes or otherwise
    case 'anycomb'
        iter = 100;
    case 'pairs'
        iter = 1000;
    otherwise
        error('Enter either ''anycomb'' or ''pairs'' for compopt')
end

stimtype = exp_site_nedata.stim;
stimtype = regexp(stimtype, 'rn\d{1,2}','match','once');
df = exp_site_nedata.df;

curdrive = gcdr;

folder = 'Ripple_Noise';
subfolder = 'downsampled_for_MID';

if isfield(exp_site_nedata, 'stimlength')
    stimlength = exp_site_nedata.stimlength;
    if df <=10
        stimmatfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_*matrix.mat',stimtype,stimlength,df)),1);
    else
        stimmatfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt10_*matrix.mat',stimtype,stimlength)),1);
    end
else

    if df <= 10
        stimmatfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-10min_DFt%d_*matrix.mat',stimtype,df)),1);
    else
        stimmatfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-10min_DFt10_*matrix.mat',stimtype)),1);
    end
end

load(stimmatfile{1});
  
% initialize correlation structure array with 10x length data (excess rows
% will be removed later)
corrstruct(length(spkstruct)*10).neuron = [];
c = 1;

rciter = 20; % real correlation iteration

% get real correlation values
for i = 1:length(spkstruct)
    
    fprintf('\nCalculating STA similarity for neuron %d...\n',spkstruct(i).neuron)
    
%     if isfield(spkstruct, 'NEs_after_thresh')
%         comb = nchoosek(spkstruct(i).NEs_after_thresh, 2);
%     else
        comb = nchoosek(spkstruct(i).NEs, 2);
%     end
    
    for j = 1:size(comb,1)
        
        fprintf('\nCalculating combination %d of %d...\n', j, size(comb,1))
        
%         if isfield(spkstruct, 'NEs_after_thresh')
%             combidx = arrayfun(@(x) find(x == spkstruct(i).NEs_after_thresh), comb(j,:));
%         else
            combidx = arrayfun(@(x) find(x == spkstruct(i).NEs), comb(j,:));
%         end
        corrstruct(c).neuron = spkstruct(i).neuron;
        corrstruct(c).NEs = comb(j,:);
        corrstruct(c).subspktrain = spkstruct(i).NEexclusive(combidx,:);
        corrstruct(c).subspktrain_count = spkstruct(i).NEexclusive_count(combidx);
        [corrstruct(c).min_spkcount, minidx] = min(corrstruct(c).subspktrain_count);
        
        % get index of spktrain with more spkcounts
        if minidx == 1
            maxidx = 2;
        else
            maxidx = 1;
        end
        
        toremove1 = corrstruct(c).subspktrain_count(maxidx) - corrstruct(c).min_spkcount;
        minSTA = quick_calc_sta(stim_mat, corrstruct(c).subspktrain(minidx,:), nlags);
        
        subsampspktrains = zeros(rciter, size(nedata.sta_spktrain,2));
        
        % sub sample group with more spikes 
        for k = 1:rciter            
            subsampspktrains(k,:) = sub_sample_spktrain(corrstruct(c).subspktrain(maxidx,:), toremove1);
        end
        
        maxSTA = quick_calc_sta(stim_mat, subsampspktrains, nlags);
        
        switch staopt
            case 'norm'
                minSTA = ne_sig_sta_from_stim_obs_resp(minSTA, corrstruct(c).subspktrain(minidx,:), stim_mat, 50, nlags, 90);
                maxSTA = ne_sig_sta_from_stim_obs_resp(maxSTA, subsampspktrains, stim_mat, 50, nlags, 90);
        end
        
        corrstruct(c).real_corrvals = corr(minSTA', maxSTA')';
        
        excspktrain = sum(corrstruct(c).subspktrain,1); %NE-exclusive spikes only
        count = sum(excspktrain);
        
        switch compopt
            
            case 'pairs'                

%                 toremove2 = count - 2*min(corrstruct(c).min_spkcount);
                train1 = zeros(iter,length(excspktrain));
                train2 = zeros(iter,length(excspktrain));

                for k = 1:iter
                    temptrain = sub_sample_spktrain(corrstruct(c).subspktrain(maxidx,:), toremove1);
%                     temptrain = sub_sample_spktrain(excspktrain, toremove2);
                    temptrain = corrstruct(c).subspktrain(minidx,:) + temptrain;  
                    splittrain = split_spktrain(temptrain);
                    train1(k,:) = splittrain(1,:);
                    train2(k,:) = splittrain(2,:);
                end
                
                stamat1 = quick_calc_sta(stim_mat, train1, nlags);
                stamat2 = quick_calc_sta(stim_mat, train2, nlags);
                
                switch staopt
                    case 'raw'
                        stasig1 = stamat1;
                        stasig2 = stamat2;
                    case 'norm'
                        [stasig1, ~] = ca_sig_sta_from_stim_obs_resp(stamat1, train1(:,nedata.nlags:end), stim, 50, 90);
                        [stasig2, ~] = ca_sig_sta_from_stim_obs_resp(stamat2, train2(:,nedata.nlags:end), stim, 50, 90);
                end
                
                corrstruct(c).surr_corrvals = diag(corr(stasig1',stasig2'));
            
            case 'anycomb'
                
                toremove2 = count - corrstruct(c).min_spkcount;
                sub_sample = zeros(iter,length(excspktrain));

                for k = 1:iter
                    sub_sample(k,:) = sub_sample_spktrain(excspktrain, toremove2);
                end
                
                sub_sample = sub_sample(:,nedata.nlags:end);
                sub_stamat = sub_sample * stim;

                switch staopt
                    case 'raw'
                        substa_sig = sub_stamat;
                    case 'norm'
                        [substa_sig, ~] = ca_sig_sta_from_stim_obs_resp(sub_stamat, sub_sample, stim, 50, 90);
                end

                sub_corrvals = triu(corr(substa_sig'),1);
                corrstruct(c).surr_corrvals = sub_corrvals(sub_corrvals ~=0);
        end
        
        c = c+1;
        
    end
end

fprintf('\n')

corrstruct(c:end) = []; %remove empty rows


