function varargout = ne_plot_ensemble_members_sta_from_sigsta(sig_sta, sigopt, normopt, saveopt, x0)

% Plots the stem plot of each independent component with the threshold for
% membership, each cNE's STA, and all member neurons' STAs.
%
% sig_sta: significant_sta struct array. See ne_get_sig_sta_structarray.m.
% sigopt: To plot normalized STAs (only significant pixels), enter True. To
% plot raw STAs, enter False.
% normopt: To normalize colors across all cNEs, enter 'together'.
% Otherwise, enter 'separate'.
% saveopt: To save figures in PDFs, enter True.
% x0: Draws dashed lines at x0 and x0+25 for choosing of x0 for MID.
% Default is [], i.e. no lines drawn.
% 
% Written 6/1/18 by JS.

if ~exist('x0','var')
    x0 = [];
end

if ~exist('saveopt', 'var')
    saveopt = 0;
end

if ~exist('normopt','var')
    normopt = 'separate';
end


load('I:\Ripple_Noise\downsampled_for_MID\rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt1_DFf5_param.mat', 'faxis')
ytick = 20:20:length(faxis);
ylab = round(faxis(20:20:length(faxis))/1000);

if sigopt
    neuron_sta = {sig_sta.neuron_sig_sta};
    NE_sta = {sig_sta.NE_sig_sta};
else
    neuron_sta = {sig_sta.neuron_sta};
    NE_sta = {sig_sta.NE_sta};
end

if saveopt
    pdffiles = cell(length(sig_sta),1);
end

for i = 1:length(sig_sta)
    
    if isempty(x0)
        
        ne_plot_ensemble_members_sta(sig_sta(i).IC_weights, sig_sta(i).CI,...
            sig_sta(i).neurons, neuron_sta{i}, NE_sta{i}, sig_sta(i).NE, normopt,...
            'NEcount', sig_sta(i).NE_event_count, 'neuroncount',...
            sig_sta(i).neuron_spike_count, 'ytick', ytick, 'ylab', ylab,...
            'filename', sig_sta(i).filename);
        
    else
        
        ne_plot_ensemble_members_sta(sig_sta(i).IC_weights, sig_sta(i).CI,...
            sig_sta(i).neurons, neuron_sta{i}, NE_sta{i}, sig_sta(i).NE, normopt,...
            'NEcount', sig_sta(i).NE_event_count, 'neuroncount',...
            sig_sta(i).neuron_spike_count, 'ytick', ytick, 'ylab', ylab,...
            'filename', sig_sta(i).filename, 'x0', x0(i));
    end
        
    
    if saveopt
        
        pdffiles{i} = sprintf('temp%d.pdf', i);
        
        export_fig(gcf, pdffiles{i}, '-nocrop');
        close all;
        
    end
        
    
end

if saveopt 
    varargout{1} = pdffiles;
end