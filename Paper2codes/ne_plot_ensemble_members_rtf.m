function ne_plot_ensemble_members_rtf(exp_site_nedata)

% Plots the RTF of each cell assembly with the RTFs of its member neurons.

% Written 5/30/18 by JS

nedata = exp_site_nedata.nedata;
NEmembers = nedata.NEmembers;
ICwt = nedata.Patterns;
neuron_rtfmat = nedata.neuron_rtfmat;
NE_rtfmat = nedata.NE_rtfmat;
CI = nedata.CI;
tmf = nedata.tmf;
xmf = nedata.xmf;

cmap = brewermap(1000,'OrRd');


for i = 1:size(NE_rtfmat,1)

    figure('units','normalized','outerposition',[0 0 1 1]);

    num_mem = length(NEmembers{i});
    
    if num_mem == 1
        nrows = 3;
        ncols = 1;
    elseif num_mem == 2
        nrows = 3;
        ncols = 2;
    elseif num_mem <= 4
        nrows = 4;
        ncols = 2;
    elseif num_mem <= 6
        nrows = 4;
        ncols = 3;
    elseif num_mem <= 9
        nrows = 5;
        ncols = 3;
    elseif num_mem <= 12
        nrows = 5;
        ncols = 4;
    elseif num_mem <= 16
        nrows = 6;
        ncols = 4;
    elseif num_mem <= 20
        nrows = 7;
        ncols = 5;
    elseif num_mem <= 25
        nrows = 8;
        ncols = 5;
    else
        error('Too many members right now.');
    end
      
       
    subplot(nrows, ncols, 1:ncols);
    stem(ICwt(:,i));
    xlim([0 size(ICwt,1)+1]);
    line([0 size(ICwt,1)+1],[CI(1) CI(1)],'Color','r','LineStyle','--');
    line([0 size(ICwt,1)+1],[CI(2) CI(2)],'Color','r','LineStyle','--');
    title(sprintf('Weights of ICs, NE #%d',i));
    xlabel('Neuron')
    ylabel('IC weight')
    
    axNE = subplot(nrows, ncols, ncols + ceil(ncols/2));
    rfmat = reshape(NE_rtfmat(i,:), length(xmf), length(tmf));

    imagesc(tmf, xmf, rfmat);
    set(axNE,'ydir', 'normal');
    tickpref;
    colormap(axNE, cmap);
    
    title(sprintf('NE #%d',i));
    xlabel('Temporal Modulation (Hz)')
    ylabel('Spectral Modulation (cyc/oct)')
    
    ax = zeros(length(NEmembers{i}),1);
    
    for j = 1:length(NEmembers{i})
        ax(j) = subplot(nrows, ncols, 2*ncols+j);
        rfmat = reshape(neuron_rtfmat(NEmembers{i}(j),:), length(xmf), length(tmf));      

        imagesc(tmf, xmf, rfmat);
        set(ax(j),'ydir', 'normal');
        tickpref;
        colormap(ax(j), cmap);

        title(sprintf('Neuron #%d',NEmembers{i}(j)));
        xlabel('Temporal Modulation (Hz)')
        ylabel('Spectral Modulation (cyc/oct)')
        
    end
    
    
end

end