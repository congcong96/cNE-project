function rtfcell = ne_plot_shared_neurons_rtf(exp_site_nedata, data, sigopt, plotopt)

if nargin == 2
    sigopt = 0;
    plotopt = 1;
elseif nargin == 3
    plotopt = 1;
end

nedata = exp_site_nedata.nedata;
nf = nedata.nf;
nlags = nedata.nlags;
nonunique = [data.neuron];
% snt = params.snt;
assemidx = [data.assemidx];
spkcount = [data.spkcount];
stacell = [data.stamat];
spkcell = [data.spkcell];


stimtype = exp_site_nedata.stim;
stimtype = regexp(stimtype, 'rn\d{1,2}','match','once');
df = exp_site_nedata.df;

if df > 10
    df = 10;
end

stimlength = exp_site_nedata.stimlength;

curdrive = gcdr;

folder = 'Ripple_Noise';
subfolder = 'downsampled_for_MID';

stimparamfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-10min_DFt%d_DFf5_param.mat',stimtype,df)),1);
load(stimparamfile{1}, 'faxis', 'MaxFM', 'MaxRD');

taxis = 0.005:0.005:0.1;

if sigopt == 1

    stamat =  cell2mat(stacell);
    spkmat = cell2mat(spkcell);
%     spkmat = spkmat(:,nlags:end);
    
    stimmatfile =  gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_matrix.mat',stimtype,stimlength,df)), 1);

    load(stimmatfile{1});

    [sta_sig, ~] = ca_sig_sta_from_stim_obs_resp(stamat, spkmat, stim_mat, 20, nlags, 95);
    cellsize = cellfun(@(x) size(x,1),spkcell);
    stacell = mat2cell(sta_sig,cellsize,size(stamat,2));
end

rtfcell = cell(length(nonunique),1);

for i = 1:length(stacell)
    
    if plotopt == 1
    
        stafig = figure;
        rtffig = figure;

        ncells = size(stacell{i},1);

        if ncells == 2
            nrows = 1;
            ncols = 2;
        elseif ncells <= 4
            nrows = 2;
            ncols = 2;
        elseif ncells <= 6
            nrows = 3;
            ncols = 2;
        elseif ncells <= 12
            nrows = 3;
            ncols = 4;
        else % ncells <= 16
            nrows = 4;
            ncols = 4;
        end
    end
    rtfcell{i} = zeros(size(stacell{i},1),16*31);

    for j = 1:size(stacell{i},1)
        

        rfmat1 = reshape(stacell{i}(j,:), nf, nlags);
               
        [tmf, xmf, rtf] = sta2rtf(rfmat1, taxis, faxis, MaxFM, MaxRD);
        
        rtfcell{i}(j,:) = (rtf(:))';
        
%         rtf_sig(j,  rtf(j,:) < rtfsiglevel (j,1) | rtf(j,:) > rtfsiglevel(j,2)) = ...
%             rtf(j, rtf(j,:) < rtfsiglevel (j,1) | rtf(j,:) > rtfsiglevel(j,2));
        
        if plotopt == 1            
            
            ytick = 10:10:length(faxis);
            ylab = round(faxis(10:10:length(faxis))/1000);
            
            figure(stafig);
            subplot(nrows,ncols,j);
            imagesc(fliplr(rfmat1));
            
            minmin = min(min(rfmat1));
            maxmax = max(max(rfmat1));
            boundary = max([abs(minmin) abs(maxmax)]);
            set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
            
            xlabel('time before spike (ms)')
            ylabel('frequency (kHz)')
            set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
            set(gca, 'xdir','reverse')
            set(gca,'ydir', 'normal');
            set(gca, 'ytick',ytick,'yticklabel',ylab)
            set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
            cmap = cschemes('rdbu', 21);
            colormap(cmap);            
            title(sprintf('Neuron #%d ONLY with cNE #%d\n(number of spikes: %d)', ...
                nonunique(i), assemidx{i}(j), spkcount{i}(j)));
            
            
            figure(rtffig);
            subplot(nrows, ncols, j);
            imagesc(tmf,xmf,rtf);
%             
%             minmin = min(min(rtf));
%             maxmax = max(max(rtf));
%             boundary = max([abs(minmin) abs(maxmax)]);
%             set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
            
            axis xy
            set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
            colormap jet;

            xlabel('temporal modulation (Hz)')
            ylabel('spectral modulation (cyc/oct)')            
            
            title(sprintf('Neuron #%d ONLY with cNE #%d\n(number of spikes: %d)', ...
                nonunique(i), assemidx{i}(j), spkcount{i}(j)));
            
        end

    end
    
    stafig.Position = [100 100 1250 500];
    rtffig.Position = [100 100 1250 500];    
    
end