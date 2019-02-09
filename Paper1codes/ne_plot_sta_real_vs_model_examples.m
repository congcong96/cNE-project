function ne_plot_sta_real_vs_model_examples(realsta, staPred, stim, neuronnum)

if nargin == 2
    neuronnum = [];
end

paramfolder = 'I:\Ripple_Noise\downsampled';
paramfile = sprintf('%s-*_DFt%d_DFf5_param.mat',stim,10);
paramfilename = gfn(fullfile(paramfolder,paramfile));
load(paramfilename{1}, 'faxis')
ytick = 20:20:length(faxis);
ylab = round(faxis(20:20:length(faxis))/1000);


fn = fieldnames(staPred);

if ~isempty(neuronnum)
    %get min and max
    realvec = realsta(neuronnum,:);
    predmat = structfun(@(x) x(neuronnum,:), staPred, 'UniformOutput', 0);
    predmat = cell2mat(struct2cell(predmat));
%     cmin = min(min([realvec;predmat]));
%     cmax = max(max([realvec;predmat]));
    
    figure;
    subplot(2,2,1)
    plot_strf_symmetric_colormap(reshape(realvec, 64, 20));
    minmin = min(realvec);
    maxmax = max(realvec);
    boundary = max(abs([minmin maxmax]));
    xlabel('Time before spike (ms)')
    ylabel('Frequency (kHz)')
    set(gca, 'xtick',0:5:20, 'xticklabel',fliplr(0:25:100))
    set(gca,'ydir', 'normal');
    set(gca, 'ytick',ytick,'yticklabel',ylab)
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);

    
    for i = 1:length(fn)
        subplot(2,2,i+1)
        plot_strf_symmetric_colormap(reshape(predmat(i,:), 64, 20));
        minmin = min(predmat(i,:));
        maxmax = max(predmat(i,:));
        boundary = max(abs([minmin maxmax]));
        xlabel('Time before spike (ms)')
        ylabel('Frequency (kHz)')
        set(gca, 'xtick',0:5:20, 'xticklabel',fliplr(0:25:100))
        set(gca,'ydir', 'normal');
        set(gca, 'ytick',ytick,'yticklabel',ylab)
        set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
        set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    end
    
end


