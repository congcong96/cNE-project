function ne_plot_ensemble_members_sta(ICwt, CI, NEmembers, neuronsta, NEsta, NEnum, varargin)

% Plots the STA of each cell assembly with the STAs of its member neurons 
% as determined via the permutation tests in ca_calc_ICA_threshold. 
%
%   nedata:     needs to have nedata.assembly_members. To obtain that, get
%               threshold from ca_calc_ICA_threshold.m first, and then run
%               ca_identify_assembly_members.m.
%   
%   nf:         #frequencies in the STA
% 
%   nlags:      #time bins in the STA 

% Updated 12/7/17 by JS to use new method of cNE STA calculation
% Updated 6/1/18 by JS to generalize it and use for different struct arrays

ip = inputParser;
addOptional(ip, 'normopt', 'separate', @(x) strcmp(x, 'together') || strcmp(x, 'separate'))
addParameter(ip, 'NEcount', [], @isscalar);
addParameter(ip, 'neuroncount', [], @(x) isvector(x) || isscalar(x));
addParameter(ip, 'ytick', [], @isvector);
addParameter(ip, 'ylab', [], @isvector);
addParameter(ip, 'nf', 64, @isscalar);
addParameter(ip, 'nlags', 20, @isscalar);
addParameter(ip, 'filename', [], @ischar);
addParameter(ip, 'x0', [], @(x) isscalar(x) || isempty(x))
parse(ip, varargin{:});

normopt = ip.Results.normopt;
NEcount = ip.Results.NEcount;
neuroncount = ip.Results.neuroncount;
ytick = ip.Results.ytick;
ylab = ip.Results.ylab;
nf = ip.Results.nf;
nlags = ip.Results.nlags;
filename = ip.Results.filename;
x0 = ip.Results.x0;
    
cmap = flipud(brewermap(1000,'rdbu'));

figure('Units', 'inches', 'Position', [5 1 8.5 8])

num_mem = length(NEmembers);

if num_mem <= 2
    nrows = 2;
    ncols = 2;
elseif num_mem <= 4
    nrows = 3;
    ncols = 2;
elseif num_mem <= 6
    nrows = 3;
    ncols = 3;
elseif num_mem <= 9
    nrows = 4;
    ncols = 3;
elseif num_mem <= 12
    nrows = 5;
    ncols = 3;
elseif num_mem <= 16
    nrows = 5;
    ncols = 4;
elseif num_mem <= 20
    nrows = 6;
    ncols = 4;
elseif num_mem <= 25
    nrows = 6;
    ncols = 5;
elseif num_mem <= 30
    nrows = 6;
    ncols = 6;
elseif num_mem <= 36
    nrows = 7;
    ncols = 6;
else
    error('Too many members right now.');
end


subplot(nrows, ncols, 1:ncols - 1);
stem(ICwt);
xlim([0 size(ICwt,1)+1]);
line([0 size(ICwt,1)+1],[CI(1) CI(1)],'Color','r','LineStyle','--');
line([0 size(ICwt,1)+1],[CI(2) CI(2)],'Color','r','LineStyle','--');
title(sprintf('Weights of ICs, cNE #%d',NEnum));
xlabel('Neuron')
ylabel('IC weight')
set(gca, 'tickdir', 'out', 'TickLength', [0.005 0.005]);


if strcmp(normopt, 'together')
    NEsta = NEsta ./ NEcount;
    neuronsta = neuronsta ./ neuroncount(:);
    minmin = min(min([NEsta; neuronsta]));
    maxmax = max(max([NEsta; neuronsta]));
    boundary = max([abs(minmin) abs(maxmax)]);    
end

axNE = subplot(nrows, ncols, ncols);
rfmat = reshape(NEsta, nf, nlags);

if strcmp(normopt, 'separate')
    minmin = min(min(rfmat));
    maxmax = max(max(rfmat));
    boundary = max([abs(minmin) abs(maxmax)]);
end

imagesc(rfmat);

if ~isempty(x0)
    line([1 nlags],[x0-0.5 x0-0.5],'Color','r','LineStyle','--');
    line([1 nlags],[x0+25.5 x0+25.5],'Color','r','LineStyle','--');
end

% xlim([0.5, nlags+0.5])
% ylim([0.5, nf+0.5])

set(axNE,'ydir', 'normal');
colormap(axNE, cmap);
set(axNE, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);

title(sprintf('cNE #%d',NEnum));
% xlabel('Time before spike (ms)')
% ylabel('Frequency (kHz)')
set(gca, 'xtick',0:5:20, 'xticklabel',[]) %fliplr(0:25:100))
set(gca,'ydir', 'normal');
if ~isempty(ytick) && ~isempty(ylab)
    set(gca, 'ytick', ytick,'yticklabel', ylab)
end
tickpref;

ax = zeros(length(NEmembers),1);

for j = 1:length(NEmembers)
    ax(j) = subplot(nrows, ncols, ncols+j);
    if size(neuronsta, 1) == length(NEmembers)
        rfmat = reshape(neuronsta(j,:), nf, nlags);
    else
        rfmat = reshape(neuronsta(NEmembers(j),:), nf, nlags);      
    end
    
    if strcmp(normopt, 'separate')
        minmin = min(min(rfmat));
        maxmax = max(max(rfmat));
        boundary = max([abs(minmin) abs(maxmax)]);
    end

    imagesc(rfmat);
    
    if ~isempty(x0)
        line([1 nlags],[x0-0.5 x0-0.5],'Color','r','LineStyle','--');
        line([1 nlags],[x0+25.5 x0+25.5],'Color','r','LineStyle','--');
    end
%     
%     xlim([0.5, nlags+0.5])
%     ylim([0.5, nf+0.5])
    
    set(ax(j),'ydir', 'normal');
    colormap(ax(j), cmap);
    set(ax(j), 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);

    title(sprintf('Neuron #%d',NEmembers(j)));
    if j > length(NEmembers) - ncols
        xlabel('Time before spike (ms)')
        set(ax(j), 'xtick',0:5:20, 'xticklabel',fliplr(0:25:100))
    else
        set(ax(j), 'xtick',0:5:20, 'xticklabel', [])
    end
    
    if mod(ncols+j, ncols) == 1
        ylabel('Frequency (kHz)')
        if ~isempty(ytick) && ~isempty(ylab)
            set(gca, 'ytick', ytick,'yticklabel', ylab)
        end
    else
        if ~isempty(ytick)
            set(ax(j), 'ytick', ytick, 'yticklabel', [])
        end
    end
    
    set(ax(j),'ydir', 'normal');
    
    tickpref;
%     if ~isempty(ylimopt)
%         ylim(ylimopt)
%     end
end

print_mfilename(mfilename);

if ~isempty(filename)
    t = suptitle(sprintf('%s-cNE%d', filename, NEnum));
    set(t, 'Interpreter', 'none')
end
    
% end % (for i)

end

