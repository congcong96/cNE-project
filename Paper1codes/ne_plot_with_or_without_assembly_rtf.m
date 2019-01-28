function rtfstruct = ne_plot_with_or_without_assembly_rtf(exp_site_nedata, plotopt)

% Updated 4/6/18 by JS to calculate RTF without plotting and clean up code

if ~exist('plotopt','var')
    plotopt = 0;
end

[stamat,totalspk,~] = ne_plot_with_or_without_assembly_sta(exp_site_nedata);

stimtype = exp_site_nedata.stim;
stimtype = regexp(stimtype, 'rn\d{1,2}','match','once');
if exp_site_nedata.df <= 10
    df = exp_site_nedata.df;
else
    df = 10;
end

if isfield(exp_site_nedata, 'stimlength')
    stimlength = exp_site_nedata.stimlength;
else
    stimlength = 10;
end

curdrive = gcdr;

folder = 'Ripple_Noise';
subfolder = 'downsampled_for_MID';

stimmatfile =  gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_matrix.mat',stimtype,stimlength,df)),1);
load(stimmatfile);
stimparamfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_param.mat',stimtype,stimlength,df)),1);
load(stimparamfile, 'faxis');

taxis = 0.005:0.005:0.1;
MaxFM = 40;
MaxRD = 4;

nf = exp_site_nedata.nedata.nf;
nlags = exp_site_nedata.nedata.nlags;
nedata = exp_site_nedata.nedata;

NEmembers = nedata.NEmembers;
numneurons = sum(cellfun('length',NEmembers));
rtfstruct(numneurons).exp = [];

rtf1 = cell(length(stamat),1);
rtf2 = cell(length(stamat),1);

c = 1;

for i = 1:length(stamat)
    
    for j = 1:size(stamat{i},1)/2
        
        rtfstruct(c).exp = exp_site_nedata.exp;
        rtfstruct(c).site = exp_site_nedata.site;
        rtfstruct(c).stim = exp_site_nedata.stim;
        rtfstruct(c).stimlength = exp_site_nedata.stimlength;
        rtfstruct(c).neuron = NEmembers{i}(j);
        rtfstruct(c).NE = i;
        
        rfmat1 = reshape(stamat{i}(j*2-1,:), nf, nlags);        
        [tmf, smf, rtf1{i}{j}] = sta2rtf(rfmat1, taxis, faxis, MaxFM, MaxRD);
        [~, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf1{i}{j}, tmf, smf);
        if any(isnan(tmtf)) || any(isnan(smtf))
            continue
        end
        [rtfstruct(c).tbmf_wo, rtfstruct(c).tbw6db_wo, rtfstruct(c).tbw3db_wo] = mtf_bmf_bw_wc(tmf_mtf, tmtf);
        [rtfstruct(c).sbmf_wo, rtfstruct(c).sbw6db_wo, rtfstruct(c).sbw3db_wo] = mtf_bmf_bw_wc(smf_mtf, smtf);
        
        
        rfmat2 = reshape(stamat{i}(j*2,:), nf, nlags);
        [tmf, smf, rtf2{i}{j}] = sta2rtf(rfmat2, taxis, faxis, MaxFM, MaxRD);
        [~, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf2{i}{j}, tmf, smf);
        [rtfstruct(c).tbmf_with, rtfstruct(c).tbw6db_with, rtfstruct(c).tbw3db_with] = mtf_bmf_bw_wc(tmf_mtf, tmtf);
        [rtfstruct(c).sbmf_with, rtfstruct(c).sbw6db_with, rtfstruct(c).sbw3db_with] = mtf_bmf_bw_wc(smf_mtf, smtf);
        
        c = c+1;
    end
end

if plotopt
    
    for i = 1:length(stamat)

        figure;

        ncells = size(stamat{i},1);

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

        for j = 1:size(stamat{i},1)/2

            n1 = 2*j-1;


            if j > 8 && mod(j,8)==1
                figure;
            end

            if n1 > 16
                n1 = mod(n1,16);
            end

            n2 = n1 + 1;

    %         rfmat1 = reshape(stamat{i}(j*2-1,:), nf, nlags);
    %         rfmat2 = reshape(stamat{i}(j*2,:), nf, nlags);
    %         
    %         [tmf, xmf, rtf1] = sta2rtf(rfmat1, taxis, faxis, MaxFM, MaxRD);
    %         [~,~,rtf2] = sta2rtf(rfmat2, taxis, faxis, MaxFM, MaxRD);

    %         minmin1 = min(min(rtf1));
    %         maxmax1 = max(max(rtf1));
    % 
    %         minmin2 = min(min(rtf2));
    %         maxmax2 = max(max(rtf2));
    %         
    %         minmin = min([minmin1 minmin2]);
    %         maxmax = max([maxmax1 maxmax2]);



            subplot(nrows, ncols, n1);
            imagesc(tmf,smf,rtf1{i}{j});
            axis xy
    %         xlabel('time before spike (ms)')
    %         ylabel('frequency (kHz)'

            set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    %         set(gca, 'clim', [0.95*minmin-eps 1.05*maxmax+eps]);
            colormap jet;

            title(sprintf('Neuron #%d without assembly\n(number of spikes: %d)',...
                nedata.NEmembers{i}(j),totalspk{i}(j*2-1)));
    %             title(sprintf('Neuron #%d non-CA',...
    %                 cadata.assembly_members{i}(j)), 'FontSize',16);

            subplot(nrows,ncols, n2);
            imagesc(tmf,smf,rtf2{i}{j});
            axis xy
    %         xlabel('time before spike (ms)')
    %         ylabel('frequency (kHz)')

            set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    %         set(gca, 'clim', [0.95*minmin-eps 1.05*maxmax+eps]);
            colormap jet;
            title(sprintf('Neuron #%d with assembly #%d\n(number of spikes: %d)', ...
                nedata.NEmembers{i}(j), i,...
                totalspk{i}(j*2)));
    %             title(sprintf('Neuron #%d CA', ...
    %                 cadata.assembly_members{i}(j)),'FontSize',16);
        end
    end
end
        

        
        
        
        
        