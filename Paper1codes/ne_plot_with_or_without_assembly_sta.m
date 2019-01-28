function [stamat,totalspk,direct] = ne_plot_with_or_without_assembly_sta(exp_site_nedata, varargin)

% Plots STA of members of the assembly when they fire and when they do not
% fire together with the assembly.

%   exp_site_cadata: standard input from cell assembly calculations.

%   opt: takes value of 0 or 1. If opt is 1, the STAs are 'normalized' 
%   between each pair of STAs per neuron by randomly removing spikes from 
%   the condition with more spikes to match the condition with less spikes.
%   If opt is 1, colors of STA plots are also normalized between each pair, i.e.
%   the colors represent the absolute values between each pair. If opt is
%   0, the STAs are not normalized and the colors are not normalized.
%   Default value is 0 (not normalized).

%   Jermyn See, updated 4/5/18.


ip = inputParser;
addRequired(ip,'exp_site_nedata', @isstruct);
addParameter(ip, 'normopt', 0, @(x) x == 0 || x == 1)
addParameter(ip, 'sigopt', 0, @(x) x == 0 || x == 1)
addParameter(ip, 'plotopt', 0, @(x) x == 0 || x == 1)
addParameter(ip, 'NEopt', [], @(x) all(x > 0 & mod(x, 1) == 0))
parse(ip, exp_site_nedata, varargin{:})

exp_site_nedata = ip.Results.exp_site_nedata;
normopt = ip.Results.normopt;
sigopt = ip.Results.sigopt;
plotopt = ip.Results.plotopt;
NEopt = ip.Results.NEopt;
   

%get cell assembly spk logicals
nedata = exp_site_nedata.nedata;
[assem_train] = ne_get_neuronal_ensemble_spktrain(exp_site_nedata, ...
    'memneuopt', 'none', 'threshalpha', 99.5, 'method', 'repeat');

NEmembers = nedata.NEmembers;

if ~isempty(NEopt)
    assem_train = assem_train(NEopt, :);
    NEmembers = NEmembers(NEopt);
end

%initialize cells
spkcell = cell(size(assem_train,1),1);
totalspk = cell(size(assem_train,1),1);
direct = cell(size(assem_train,1),1);

for i = 1:size(assem_train,1)

    direct{i} = zeros(length(NEmembers{i}),1);
    
    for j = 1:length(NEmembers{i})
        
        neuron = NEmembers{i}(j);
        if isfield(nedata, 'sta_spktrain')
            spktrain = nedata.sta_spktrain(neuron,:);
        else
            spktrain = nedata.spktrain(neuron,:);
        end
       
        compmatrix = [assem_train(i,:);spktrain];
        with_idx = (compmatrix(2,:) >= 1 & compmatrix(1,:) == 1);
        without_idx = (compmatrix(2,:) >= 1 & compmatrix(1,:) == 0);       
        
        
        without_spktrain = spktrain .* without_idx;
        with_spktrain = spktrain .* with_idx;
        
        %normalization starts here
        
        if normopt == 1
            wo_spkcount = sum(without_spktrain);
            w_spkcount = sum(with_spktrain);
            spkc_diff = wo_spkcount - w_spkcount;
            
            %if num with assembly spikes > num without assembly spikes
            if spkc_diff < 0
                direct{i}(j) = 0;
                with_spktrain(without_idx) = 0;
                with_spktrain = sub_sample_spktrain(with_spktrain, abs(spkc_diff));

            %if num without assembly spikes > num with assembly spikes
            elseif spkc_diff > 0
                direct{i}(j) = 1;
                without_spktrain(with_idx) = 0;
                without_spktrain = sub_sample_spktrain(without_spktrain, abs(spkc_diff));

                
            end
                        
            
        end
        
        %normalization ends here
        
        spkcell{i}(2*(j)-1,:) = without_spktrain;
        spkcell{i}(2*(j),:) = with_spktrain;
        
        
    end
    totalspk{i} = sum(spkcell{i},2);
end

spkmat = cell2mat(spkcell);

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
load(stimmatfile{1});
stimparamfile = gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_param.mat',stimtype,stimlength,df)),1);
load(stimparamfile{1}, 'faxis');
ytick = 10:10:length(faxis);
ylab = round(faxis(10:10:length(faxis))/1000);

    
[stim, resp] = ne_create_stim_trial_from_stim_matrix(stim_mat, spkmat, nedata.nlags);

stamat_sep = resp * stim;

[nf, ntrials] = size(stim_mat);
nlags = size(stamat_sep,2) / nf; 
% ca_plot_cell_assembly_stamat(stamat_sep, nf, nlags, 'STA');

if sigopt == 1
    [stamat_sep, ~] = ne_sig_sta_from_stim_obs_resp(stamat_sep, spkmat, stim_mat, 20, nlags, nf, 95);
end

count = 1;
for ii = 1:length(spkcell)
    nn = size(spkcell{ii},1);
    stamat{ii} = stamat_sep(count:count+nn-1,:);
    count = count+nn;
end

if plotopt
    
    cmap = brewmaps('rdbu', 1000);
    
    for iii = 1:length(spkcell)
        figure;
        colormap(cmap);

        ncells = size(spkcell{iii},1);

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

        if normopt == 0

            for iv = 1:size(spkcell{iii},1)
                n = iv;

                if iv > 16 && mod(iv,16)==1
                    figure;
                    colormap(cmap);
                end

                if n > 16
                    n = mod(n,16);
                    if n == 0
                        n = 16;
                    end
                end

                subplot(nrows, ncols, n);
                quick_plot_sta(reshape(stamat{iii}(iv,:),nf,nlags));
                
                boundary = max(abs(stamat{iii}(iv,:)));
                set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);

%                 rfmat = fliplr(reshape(stamat{iii}(iv,:), nf, nlags));
% 
%                 plot_strf_symmetric_colormap(rfmat);


                if mod(iv,2) == 1
                    title(sprintf('Neuron #%d without cNE\n(number of spikes: %d)',...
                        NEmembers{iii}(iv/2+0.5),totalspk{iii}(iv)));

                else
                    title(sprintf('Neuron #%d with cNE #%d\n(number of spikes: %d)', ...
                        NEmembers{iii}(iv/2), iii,...
                        totalspk{iii}(iv)));
                end

            end

        elseif normopt == 1

            for iv = 1:size(spkcell{iii},1)/2

                n1 = 2*iv-1;


                if iv > 8 && mod(iv,8)==1
                    figure;
                    colormap(cmap);
                end

                if n1 > 16
                    n1 = mod(n1,16);
                end

                n2 = n1 + 1;

                rfmat1 = reshape(stamat{iii}(iv*2-1,:), nf, nlags);
                rfmat2 = reshape(stamat{iii}(iv*2,:), nf, nlags);


                minmin1 = min(min(rfmat1));
                maxmax1 = max(max(rfmat1));

                minmin2 = min(min(rfmat2));
                maxmax2 = max(max(rfmat2));

                minmin = min([minmin1 minmin2]);
                maxmax = max([maxmax1 maxmax2]);

                boundary = max([abs(minmin) abs(maxmax)]);



                subplot(nrows, ncols, n1);
%                 imagesc(rfmat1);
%                 xlabel('time before spike (ms)')
%                 ylabel('frequency (kHz)')
%                 cmap = cschemes('rdbu', 21);
%                 colormap(cmap);
%  
%                 set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
%                 set(gca, 'xdir','reverse')
%                 set(gca,'ydir', 'normal');
%                 set(gca, 'ytick',ytick,'yticklabel',ylab)
%                 set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
                quick_plot_sta(rfmat1);
                set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
                title(sprintf('Neuron #%d without cNE\n(number of spikes: %d)',...
                    NEmembers{iii}(iv),totalspk{iii}(iv*2-1)));


                subplot(nrows,ncols, n2);
%                 imagesc(rfmat2);
%                 xlabel('time before spike (ms)')
%                 ylabel('frequency (kHz)')
%                 colormap(cmap);
%                 set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
%                 set(gca, 'xdir','reverse')
%                 set(gca,'ydir', 'normal');
%                 set(gca, 'ytick',ytick,'yticklabel',ylab)
%                 set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
                quick_plot_sta(rfmat2);
                set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
                title(sprintf('Neuron #%d with cNE #%d\n(number of spikes: %d)', ...
                    NEmembers{iii}(iv), iii, totalspk{iii}(iv*2)));
   
            end
        end
    end
        
end % (for iii)

return
