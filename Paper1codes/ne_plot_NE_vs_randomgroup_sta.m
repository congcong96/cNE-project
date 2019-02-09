function NEgroupsta = ne_plot_NE_vs_randomgroup_sta(exp_site_nedata, NEgroupinfo, ngrps, plotopt)

if ~exist('plotopt','var')
    plotopt = 0;
end

if ~exist('ngrps','var')
    ngrps = 3;
end

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;
members = nedata.NEmembers;
stimtype = exp_site_nedata.stim;
stimtype = regexp(stimtype, 'rn\d{1,2}','match','once');

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


% get cNE member spiketrain
NEtrain = cell2mat(cellfun(@(x) sum(spktrain(x,:), 1), members, 'UniformOutput', 0));
NEspkcount = sum(NEtrain,2);

randgrps = cell(length(NEgroupinfo),1);
randgrptrain = cell(length(NEgroupinfo),1);
randgrpspkcount = cell(length(NEgroupinfo), 1);

for i = 1:length(NEgroupinfo)
    
    if isempty(NEgroupinfo(i).NE_info_extrap)
        continue
    end
    
    
    grpsize = size(NEgroupinfo(i).random_groups,1);
    
    if ngrps < grpsize        
        randgrps{i} = mat2cell(NEgroupinfo(i).random_groups(randsample(grpsize,ngrps),:),...
            ones(1, ngrps), length(NEgroupinfo(i).NEmembers));
    else
        randgrps{i} = mat2cell(NEgroupinfo(i).random_groups, ones(1,grpsize),...
            length(NEgroupinfo(i).NEmembers));
    end
    randgrptrain{i} = cell2mat(cellfun(@(x) sum(spktrain(x,:), 1), randgrps{i}, 'UniformOutput', 0));
    randgrpspkcount{i} = sum(randgrptrain{i},2)';
    
end

% initialize sta struct array
NEgroupsta(length(NEgroupinfo)).NE = [];


for j = 1:length(NEgroupinfo) 
    
    if isempty(NEgroupinfo(j).NE_info_extrap)
        continue
    end 
        
    fprintf('Calculating NE #%d vs random groups STA...\n', j);

    NEgroupsta(j).NE = j;
    NEgroupsta(j).NEmembers = members{j};
    NEgroupsta(j).randgroup = randgrps{j};

    min_spikes = min([NEspkcount(j) randgrpspkcount{j}]);

    samp_NE = sub_sample_spktrain(NEtrain(j,:), NEspkcount(j) - min_spikes);
    
    samp_randgrp = zeros(ngrps, size(spktrain,2));
    sta_randgrp = cell(ngrps,1);
    
    for i = 1:length(randgrpspkcount{j})
        samp_randgrp(i,:) = sub_sample_spktrain(randgrptrain{j}(i,:), randgrpspkcount{j}(i) - min_spikes);
        sta_randgrp{i} = ca_sta_from_locator_stimulus(samp_randgrp(i,:), stim_mat, nlags);
    end
        
    NEgroupsta(j).min_spikes = min_spikes;
    NEgroupsta(j).sta_NE = ca_sta_from_locator_stimulus(samp_NE, stim_mat, nlags);
    NEgroupsta(j).sta_randgroups = sta_randgrp;      
       
end

if plotopt == 1

    for i = 1:length(NEgroupsta)
        
        if isempty(NEgroupinfo(i).NE_info_extrap)
            continue
        end
        figure;
        
        rfmats = [NEgroupsta(i).sta_NE; NEgroupsta(i).sta_randgroups];    
        minmin = min(cellfun(@(x) min(x(:)), rfmats));        
        maxmax = max(cellfun(@(x) max(x(:)), rfmats));       

        boundary = max([abs(minmin) abs(maxmax)]);
        
        for j = 1:length(rfmats)
            
            if length(rfmats) <= 4
                nrows = 2;
                ncols = 2;
            elseif length(rfmats) <= 6
                nrows = 3;
                ncols = 2;
            elseif length(rfmats) <= 9
                nrows = 3;
                ncols = 3;
            else
                error('ngrps must be <= 8')
            end
            
            subplot(nrows,ncols,j)
            imagesc(rfmats{j}./boundary);
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
            
            if j == 1
                title(sprintf('cNE #%d\n(number of spikes: %d)',...
                    NEgroupsta(i).NE, NEgroupsta(i).min_spikes));
            else
                title(sprintf('cNE #%d - random group #%d\n(number of spikes: %d)',...
                    NEgroupsta(i).NE, j-1, NEgroupsta(i).min_spikes));
            end
            
            c = colorbar;
            c.TickDirection = 'out';
            ylabel(c, 'Normalized power')
                

%         subplot(2,3,c+3);
%         imagesc(rfmat2);
%         xlabel('time before spike (ms)')
%         ylabel('frequency (kHz)')
%         colormap(cmap);
%         set(gca, 'xtick',1:5:20, 'xticklabel',fliplr(25:25:100))
%     %     set(gca, 'xdir','reverse')
%         set(gca,'ydir', 'normal');
%         set(gca, 'ytick',ytick,'yticklabel',ylab)
%         tickpref;
%         set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%         title(sprintf('Neuron #%d\n(number of spikes: %d)', ...
%             NEgroupsta(i).neuron, NEgroupsta(i).min_spikes));


%         if c == 3 && i ~= length(NEgroupsta)
%             print_mfilename(mfilename);
%             figure;
%             c = 1;
%         else
%             c = c+1;
%         end
% 
        end
    print_mfilename(mfilename);    
    end
    
end

