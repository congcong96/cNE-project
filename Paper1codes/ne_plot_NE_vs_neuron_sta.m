function NEneuronsta = ne_plot_NE_vs_neuron_sta(exp_site_nedata, plotopt)

% plotopt: cNE that will be plotted

if ~exist('plotopt','var')
    plotopt = 0;
end

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;
neact = nedata.Activities;
nethresh = nedata.NEthresh;
members = nedata.NEmembers;
stimtype = exp_site_nedata.stim;
stimtype = regexp(stimtype, 'rn\d{1,2}','match','once');

spike_count_threshold = 300;

if isfield(exp_site_nedata, 'stimlength')
    stimlen = exp_site_nedata.stimlength;
else
    stimlen = 10;
end


if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
    spktrain = nedata.spktrain;
else
    dft = 10;
    spktrain = nedata.sta_spktrain;
end


curdrive = gcdr;

folder = 'Ripple_Noise';
subfolder = 'downsampled_for_MID';

stimmatfile =  gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_matrix.mat',stimtype,stimlen,dft)),1);
load(stimmatfile{1});
stimparamfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_param.mat',stimtype,stimlen,dft)),1);
load(stimparamfile{1}, 'faxis');
ytick = 10:10:length(faxis);
ylab = round(faxis(10:10:length(faxis))/1000);


% get neurons part of at least one cNE
neunespks = cellfun(@(x) logical(spktrain(x,:)), members, 'UniformOutput', 0);
spkcount = cellfun(@(x) sum(x, 2), neunespks, 'UniformOutput', 0);

% if exp_site_nedata.df > 10
%     neactds = neact;
%     neact = zeros(length(nethresh), size(spktrain,2));
% end
% 
% NEraster = zeros(size(neact));
% 
% %get cNE events
% for i = 1:length(nethresh)
%     
%     if exp_site_nedata.df > 10     
%         
%         temptrain = interp(neactds(i,:),exp_site_nedata.df/dft);
%         if length(temptrain) == size(spktrain,2)
%             neact(i,:) = temptrain;
%         else
%             lendiff = size(neact,2) - length(temptrain);
%             neact(i, :) = temptrain(1:end+lendiff);
%         end
%         
%     end         
%     NEraster(i,:) = neact(i,:) >= nethresh(i);
% end

NEraster = ne_upsample_NEact_using_member_neuron_activity(neact, members, spktrain, nethresh);

NEeventcount = sum(NEraster,2);


% initialize info struct array
totalneurons = sum(cellfun(@length, spkcount));
NEneuronsta(totalneurons).neuron = [];

c = 1;

for j = 1:length(spkcount)    
    
    for k = 1:length(spkcount{j})
        
        fprintf('Calculating NE #%d vs neuron #%d STA...\n', j, members{j}(k));
        
        NEneuronsta(c).neuron = members{j}(k);
        NEneuronsta(c).NE = j;
        
        min_spikes = min([NEeventcount(j) spkcount{j}(k)]);
        
        if min_spikes > spike_count_threshold       

            if min_spikes == NEeventcount(j)
                samp_NE = NEraster(j,:);            
                samp_neuron = sub_sample_spktrain(neunespks{j}(k,:), spkcount{j}(k) - min_spikes);
            else
                samp_NE = sub_sample_spktrain(NEraster(j,:), NEeventcount(j) - min_spikes);
                samp_neuron = neunespks{j}(k,:);
            end

            NEneuronsta(c).min_spikes = min_spikes;
            NEneuronsta(c).sta_NE = ca_sta_from_locator_stimulus(samp_NE, stim_mat, nlags);
            NEneuronsta(c).sta_neuron = ca_sta_from_locator_stimulus(samp_neuron, stim_mat, nlags);

            c = c+1;
        else
            continue
        end
    end
end

if plotopt > 0

    cNE = plotopt;
    staplots = NEneuronsta([NEneuronsta.NE] == cNE);

    for i = 1:length(staplots)
        
        figure;

        rfmat1 = staplots(i).sta_NE;
        rfmat2 = staplots(i).sta_neuron;

        minmin1 = min(min(rfmat1));
        maxmax1 = max(max(rfmat1));

        minmin2 = min(min(rfmat2));
        maxmax2 = max(max(rfmat2));

        minmin = min([minmin1 minmin2]);
        maxmax = max([maxmax1 maxmax2]);

        boundary = max([abs(minmin) abs(maxmax)]);

        subplot(2,1,1)
        imagesc(rfmat1./boundary);
        xlabel('time before spike (ms)')
        ylabel('frequency (kHz)')
        cmap = cschemes('rdbu', 1000);
        colormap(cmap);
        set(gca, 'xtick',1:5:20, 'xticklabel',fliplr(25:25:100))
    %     set(gca, 'xdir','reverse')
        set(gca,'ydir', 'normal');
        set(gca, 'ytick',ytick,'yticklabel',ylab)
        tickpref;
        set(gca, 'clim', [-1 1]);
        title(sprintf('NE #%d\n(number of spikes: %d)',...
            staplots(i).NE, staplots(i).min_spikes));
        c = colorbar;
        c.TickDirection = 'out';

        subplot(2,1,2);
        imagesc(rfmat2./boundary);
        xlabel('time before spike (ms)')
        ylabel('frequency (kHz)')
        colormap(cmap);
        set(gca, 'xtick',1:5:20, 'xticklabel',fliplr(25:25:100))
    %     set(gca, 'xdir','reverse')
        set(gca,'ydir', 'normal');
        set(gca, 'ytick',ytick,'yticklabel',ylab)
        tickpref;
        set(gca, 'clim', [-1 1]);
        title(sprintf('Neuron #%d\n(number of spikes: %d)', ...
            staplots(i).neuron, staplots(i).min_spikes));
        c = colorbar;
        c.TickDirection = 'out';
        
        set(gcf, 'Position', [10 10 500 900]);
        print_mfilename(mfilename);


        
%         pos2 = get(subplot(2,1,2),'Position');
%         c = colorbar('Position', [pos2(1)+pos2(3)+0.02  pos2(2)  0.02  (pos2(2)+pos2(3))]);
%         c.TickDirection = 'out';
% 
%         if c == 3 && i ~= length(staplots)
%             print_mfilename(mfilename);
%             figure;
%             c = 1;
%         else
%             c = c+1;
%         end


    end
    
end

