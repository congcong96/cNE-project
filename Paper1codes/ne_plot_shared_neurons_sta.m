function [data, varargout] = ne_plot_shared_neurons_sta(exp_site_nedata, varargin)

%spkthresh: number of spikes the smallest subgroup must have to be
%considered.

%spktrainopt: 
% 'with' - takes spikes that belong to a cNE, regardless of overlap with
%   other cNEs
% 'exclusiveonly' - takes spikes that belong exclusively to a cNE
% 'all' - takes all combinations, i.e. all spikes, cNE-exclusive spikes,
%   non-cNE spikes, all-cNE spikes
% 'exclusive/total' - takes exclusively-only spikes and all spikes


%normopt: normalizes the number of spikes for each sta based on the group
%with the lowest number of spikes. can only be 1 if spktrain opt is 1.

%staopt: eliminates noise in STAs if staopt == 1

%plotopt: if plotopt == 0, STAs will not be plotted.

%stasimopt: if 1, skips STA calculation. Should only be 1 when declaring
%this code from ne_calc_shared_neurons_sta_similarity.

% Updated by JS, 7/6/17, included inputParser and gave more options for
% spiketrainopt
% Updated by JS, 9/3/17, changed code for changes in
% ne_get_shared_neuron_spktrain_subsets.m. However, only works for
% ne_calc_shared_neurons_sta_similarity.m for now. Will fix plotting soon.
% Also added stasimopt to skip the sta calculation part when using
% ne_calc_shared_neurons_sta_similarity.m since STAs are calculated in
% there.

ip = inputParser;
addRequired(ip, 'exp_site_nedata', @isstruct)
addOptional(ip, 'spkthresh', 100, @(x) mod(x,1) == 0)
addParameter(ip, 'normopt', 0, @(x) x == 0 || x == 1)
addParameter(ip, 'stanormopt', 1, @(x) x == 0 || x == 1)
addParameter(ip, 'stasigopt', 0, @(x) x == 0 || x == 1)
addParameter(ip, 'plotopt', [], @(x) isempty(x) || strcmp(x, 'exclusiveonly') || strcmp(x, 'all'))
addParameter(ip, 'stasimopt', 0, @(x) x == 0 || x == 1)
addParameter(ip, 'neuopt', [], @(x) isempty(x) || isvector(x) || isscalar(x)) %specifies whether specific neurons are looked at

parse(ip, exp_site_nedata, varargin{:});
exp_site_nedata = ip.Results.exp_site_nedata;
spkthresh = ip.Results.spkthresh;
normopt = ip.Results.normopt;
stanormopt = ip.Results.stanormopt;
stasigopt = ip.Results.stasigopt;
plotopt = ip.Results.plotopt;
stasimopt = ip.Results.stasimopt;
neuopt = ip.Results.neuopt;

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;

% get spktrain subsets based on assembly association
spkstruct = ne_get_shared_neuron_spktrain_subsets(exp_site_nedata, 'spkthresh', spkthresh, 'neuopt', neuopt);

if isempty(spkstruct)
    data = [];
    warning('No shared neurons in this dataset')
    return
end
 
% if ~strcmp(plotopt, 'exclusiveonly')
%     normopt = 0;
% end


% % if spkthresh is a variable, use it to consider only groups with more than
% % the number of spikes specified by spkthresh
% if ~isempty(spkthresh) == 1 && strcmp(spktrainopt, 'exclusiveonly')
%     for i = 1:length(spkstruct)
%         
%         spkcount = spkstruct(i).NEexclusive_count;
%         keepidx = spkcount >= spkthresh;
%         spkstruct(i).NEs = spkstruct(i).NEs(keepidx);
%         spkstruct(i).NEexclusive = spkstruct(i).NEexclusive(keepidx,:);
%         spkstruct(i).NEexclusive_count = spkstruct(i).NEexclusive_count(keepidx);
% 
%     end
%     spkstruct(cellfun(@(x) length(x) < 2,{spkstruct.NEs})) = [];
% end

% if isempty(spkstruct)
%     data = [];
%     warning('No shared neurons in this dataset after applying threshold')
%     return
% end

if stasimopt == 1
    data = spkstruct;
    return
end

% if normopt is 1 and spktrainopt is 1, normalize with-exclusive spikes by
% the assembly with the lowest number of spikes.
% 
% if strcmp(spktrainopt, 'exclusiveonly') && normopt == 1
%     
%     count = 1;
%     niter = 100;
%     rownum = 1;
%     %get number of spiketrains being subsampled
%     spkcount = {spkstruct.NEexclusive_count};
%     s = length(cell2mat(spkcount')) - length(spkcount)*2;
%     
%     sampcell = cell(s,1);
%     tracker = zeros(s,1);
%     
%     for j = 1:length(spkstruct)        
%             
%         minspk = min(spkcount{j});
%         
%         for k = 1:length(spkcount{j})
%             
%             spkdiff = spkcount{j}(k) - minspk;             
%             if k ~= 1 && spkdiff ~= 0
%                 spkcount{j}(k) = minspk;
%                 for L = 1:niter
%                     sampcell{count}(L,:) = sub_sample_spktrain(...
%                         spkstruct{j}(k,:), spkdiff);
%                     tracker(count) = rownum;
%                 end
%                 count = count + 1;
%             end
%             rownum = rownum + 1;
% 
%         end
%     end
%    sampcell = sampcell(~cellfun('isempty',sampcell));
%    tracker = tracker(1:length(sampcell));
%        
% end

% spkcount = cellfun(@(x) sum(x,2), spkcell, 'UniformOutput', 0);

if strcmp(plotopt, 'exclusiveonly') && normopt == 1
    
    NEexctrain = {spkstruct.NEexclusive}';
    NEexccount = {spkstruct.NEexclusive_count}';
    
    spkcell = cell(length(NEexccount),1);
    titles = cell(length(NEexccount), 1);
    
    for i = 1:length(NEexccount)
        
        min_spikes = min(NEexccount{i});        
        spkcell{i} = zeros(size(NEexctrain{i}));
        titles{i} = cell(length(NEexccount{i}), 1);
        
        for j = 1:length(NEexccount{i})
            
            spkcell{i}(j,:) = sub_sample_spktrain(NEexctrain{i}(j,:), NEexccount{i}(j) - min_spikes);
            titles{i}{j} = sprintf('cNE%d-exclusive spikes, n = %d', spkstruct(i).NEs(j), spkstruct(i).NEexclusive_count(j));
            
        end
    end
    
    cellidx = cellfun(@(x) size(x, 1), spkcell);
    
elseif strcmp(plotopt, 'exclusiveonly') && normopt == 0
    
    spkcell = {spkstruct.NEexclusive}';
    cellidx = cellfun(@(x) size(x,1), spkcell);
    titles = cellfun(@(y, z) arrayfun(@(x) sprintf('cNE%d-exclusive spikes', x), y, 'UniformOutput', 0), {spkstruct.NEs}', 'UniformOutput', 0);
    
elseif strcmp(plotopt, 'all')
    
    fn = {'spktrain'; 'withoutNE'; 'withanyNE'; 'withallNEs'; 'withNE'; 'NEexclusive'};
    allfn = fieldnames(spkstruct);
    [~, idx] = setdiff(allfn, fn);
    tempstruct = rmfield(spkstruct, allfn(idx));
    tempcell = squeeze(struct2cell(tempstruct));
    
    spkcell = cell(size(tempcell,2), 1);
    stdtitles = {'All spikes'; 'cNE-i spikes'; 'cNE-a spikes'; 'cNE-all spikes'};
    withtitles = cell(size(tempcell,2), 1);
    exctitles = cell(size(tempcell,2), 1);
    
    for i = 1:size(tempcell, 2)
        spkcell{i} = cell2mat(tempcell(:,i));
        withtitles{i} = arrayfun(@(x) sprintf('cNE%d spikes', x), spkstruct(i).NEs, 'UniformOutput', 0);
        exctitles{i} = arrayfun(@(x) sprintf('cNE%d-exclusive spikes', x), spkstruct(i).NEs, 'UniformOutput', 0);
    end
    
    titles = cellfun(@(x,y) [stdtitles;x;y], withtitles, exctitles, 'UniformOutput', 0);
    cellidx = cellfun(@(x) size(x,1), spkcell);
    
end

%convert to matrix for STA calculation
spkmat = cell2mat(spkcell);

% get stimulus matrix
stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);

% get STAs
stamat = quick_calc_sta(stimstr.stimulus, spkmat, nlags);

if stanormopt
    stamat = stamat ./ sum(spkmat, 2);
end

if stasigopt == 1
    [sta_sig, ~] = ne_sig_sta_from_stim_obs_resp(stamat, locator, stim, 100, nlags, 90);
    stacell = mat2cell(sta_sig, cellidx, size(sta_sig, 2));
else %if stasigopt == 0
    stacell = mat2cell(stamat, cellidx, size(stamat, 2));
end

varargout{1} = stacell;

if plotopt
    
    cmap = cschemes('rdbu', 100);
    
    for i = 1:length(stacell)
    
        figure('Position', [582 344 871 351]);
        colormap(cmap);
        num = size(stacell{i}, 1);
        
        if num <= 2
            nrows = 1;
            ncols = 2;
        elseif num <= 4
            nrows = 2;
            ncols = 2;
        elseif num <= 6
            nrows = 2;
            ncols = 3;
        elseif num <= 9
            nrows = 3;
            ncols = 3;
        elseif num <= 12
            nrows = 4;
            ncols = 3;
        elseif num <= 16
            nrows = 4;
            ncols = 4;
        elseif num <= 20
            nrows = 5;
            ncols = 4;
        elseif num <= 25
            nrows = 5;
            ncols = 5;
        elseif num <= 30
            nrows = 5;
            ncols = 6;
        elseif num <= 36
            nrows = 6;
            ncols = 6;
        else
            warning('Too many members!');
            continue
        end
        
        if normopt || stanormopt
            boundary = max(abs(stacell{i}(:)));
        end
        
        for j = 1:size(stacell{i},1)
                      
            subplot(nrows, ncols, j)
            
            if ~normopt && ~stanormopt
                boundary = max(abs(stacell{i}(j,:)));
            end
                
            quick_plot_sta(reshape(stacell{i}(j,:), [], nlags));      
            
            title(titles{i}{j});
            
        end
        suptitle(sprintf('Neuron %d', spkstruct(i).neuron));
    end
end

data = spkstruct;
            
            

% 
% [nf, ~] = size(stimulus);
% nlags = size(sta_sig,2) / nf; 
% 
% % ca_plot_cell_assembly_stamat(stamat_sep, nf, nlags, 'STA');
% 
% stamat = cell(length(nonunique),1);
% count = 1;
% 
% if strcmp(spktrainopt, 'all')
% 
%     for ii = 1:length(nonunique)
%         nn = size(snt{ii},1)+2;
%         stamat{ii} = sta_sig(count:count+nn-1,:);
%         count = count+nn;
%     end
% 
% elseif strcmp(spktrainopt, 'with') || strcmp(spktrainopt, 'exclusive/total')
% 
%     for ii = 1:length(nonunique)
%         nn = size(snt{ii},1);
%         stamat{ii} = sta_sig(count:count+nn-1,:);
%         count = count + nn;
%     end
%     
% elseif strcmp(spktrainopt, 'exclusiveonly')
%     
%     for ii = 1:length(nonunique)
%         nn = size(snt{ii},1) - 1;
%         stamat{ii} = sta_sig(count:count+nn-1,:);
%         count = count + nn;
%     end
%     
% end
% 
% n = 1;
% data.stamat = stamat;
% 
% if ~isempty(plotopt)
% 
% 
%     for iii = 1:length(stamat)
% 
%         h1 = figure;
%         if size(stamat{iii},1) <= 4
%             nrows = 2;
%             ncols = 2;
%         elseif size(stamat{iii},1) <= 9
%             nrows = 3;
%             ncols = 3;
%         else
%             nrows = 3;
%             ncols = 4;
%         end
% 
%         for iv = 1:size(stamat{iii},1)
%             figure(h1)
%             subplot(nrows, ncols, iv);
% 
%             rfmat = reshape(stamat{iii}(iv,:), nf, nlags); 
%             imagesc(fliplr(rfmat));
%             xlabel('time before spike (ms)', 'FontSize', 14)
%             ylabel('frequency (kHz)', 'FontSize', 14)
%             set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100, 'FontSize', 12)
%             set(gca, 'xdir','reverse')
%             set(gca,'ydir', 'normal');
%             set(gca, 'ytick',ytick,'yticklabel',ylab, 'FontSize', 12)
%             set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
% %             colormap jet;
%             cmap = cschemes('rdbu', 21);
%             colormap(cmap);
% 
%             if normopt == 1
%                 temp = stamat{iii}(2:end,:);
%                 maxmax = max(temp(:));
%                 minmin = min(temp(:));
%                 boundary = max([abs(minmin) abs(maxmax)]);    
%                 set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%             end
% 
%             if strcmp(spktrainopt, 'all')
% 
%                 if iv == 1
%                     title(sprintf('Neuron #%d without assembly\n(number of spikes: %d)',...
%                         nonunique(iii),spkcount{iii}(1)));
%                 elseif iv == size(stamat{iii},1)
%                     title(sprintf('Neuron #%d total\n(number of spikes: %d)',...
%                         nonunique(iii),spkcount{iii}(size(stamat{iii},1))));
%                 elseif iv == size(stamat{iii},1)-1
%                     title(sprintf('Neuron #%d with all assemblies\n(number of spikes: %d)',...
%                         nonunique(iii),spkcount{iii}(size(stamat{iii},1)-1)));
%                 else
%                     title(sprintf('Neuron #%d ONLY with assembly #%d\n(number of spikes: %d)', ...
%                         nonunique(iii), assemidx{iii}(iv-1), spkcount{iii}(iv)));
%                 end
% 
%             elseif strcmp(spktrainopt, 'exclusiveonly')
%                
%                 title(sprintf('Neuron #%d NE #%d\n(number of spikes: %d)', ...
%                     nonunique(iii), assemidx{iii}(iv), spkcount{iii}(iv)));
%             
%             
%             elseif strcmp(spktrainopt, 'exclusive/total')
%                 
%                 if iv ~= size(stamat{iii},1)
%                     title(sprintf('Neuron #%d with NE #%d\n(number of spikes: %d)',...
%                         nonunique(iii), assemidx{iii}(iv), spkcount{iii}(iv)));
%                 else
%                     title(sprintf('Neuron #%d total STA\n(number of spikes: %d)',...
%                         nonunique(iii), spkcount{iii}(iv)))
%                 end
%             end
%             
%             h1.Position = [100 100 1250 750];
% 
%         end
%         n = n + 1;
% 
%     end % (for iii)
% 
% 
% end
             