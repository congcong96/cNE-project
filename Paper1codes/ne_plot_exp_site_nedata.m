function varargout = ne_plot_exp_site_nedata(exp_site_nedata, startx, labelopt)
% ca_plot_cell_assemblies_nedata Display cell assemblies, activities
% 
%     ca_plot_cell_assemblies_data(exp_site_nedata, startx)
%
%     exp_site_nedata : struct of data from cell assembly calculations for one
%     experiment and one site.
%
%     Has the form:
%
%     exp_site_nedata = 
% 
%            exp: experiment date
%           site: penetration number
%           stim: 'rn1', 'rn4', etc.
%             df: downsampling factor. Not overall downsampling factor
%         nedata: struct of cell assembly data
%
%
%     nedata : struct holding cell assembly results. 
%     Has the following fields:
%
%     nedata = 
% 
%           spktrain: spike train matrix
%              fsdvd: sampling rate of dvd audio system
%                 df: overall envelope downsampling factor
%           position: unit position, in um
%           Patterns: cell assemblies. Each column is an assembly
%         Activities: time course of assemblies. Each row is an assembly
%                 nf: # freq's in STA
%              nlags: # time bins in STA
%             stamat: STA of each neuron. Each row is a neuron.
%          ca_stamat: STA of each assembly. Each row is an assembly
%
%     startx : Optional. Beginning bin of raster and activity plots. Default = 0. 
%     labelopt: Optional. If 0 labels y axis as depth of neurons. If 1
%     labels y axis as neuron number. Default is 0.
% 

% Check input arguments
narginchk(1,3);
nargoutchk(0,2);

if ( nargin == 1 )
    startx = 0;
    labelopt = 0;
end

if ( nargin == 2 )
    if ( isempty(startx) )
        startx = 0;
    end
    labelopt = 0;
end

nedata = exp_site_nedata.nedata;
Activitymatrix = nedata.spktrain;
% 
% temp = Activitymatrix;
% for i = 1:15
%     Activitymatrix (2*i-1:2*i,:) = [temp(i,:);temp(i+15,:)];
% end

Patterns = nedata.Patterns;
Activities = nedata.Activities;
% position = [];
if isfield(nedata, 'position')
    position = nedata.position;
end

try
    dt = nedata.df / nedata.fsdvd * 1000; % spike train bin size, in ms
catch
    dt = nedata.binsize; 
end

if exist('position','var') == 1

    try
        labels = unique(position,'stable');
    catch
        position = cell2mat(position);
        labels = unique(position,'rows','stable');
    end
    numrepeats = zeros(length(labels),1);
end



if exist('labels','var') && iscell(labels)
    
    for i = 1:length(labels)
        index = strfind(position,labels{i});
        numrepeats(i) = sum(cell2mat(index)==1);
    end
    tick = cumsum(numrepeats);

        
elseif exist('labels','var')
        
    temp = labels;
    labels = cell(size(temp,1),1);
    
%     for i = 1:size(labels,1)
%         index = ismember(position,temp(i,:),'rows');
%         numrepeats(i) = sum(index);
% %         labels(i) = {num2str(temp(i,:))};
%         labels(i) = {num2str(temp(i,2))};
%     end
%         
    tick = cumsum(numrepeats);

end

cmap1 = cschemes('rdbu', 21);
cmap2 = cschemes('reds', 50);
cmap2 = flipud([cmap2;1,1,1]);

% Estimate pairwise channel correlations
correlationmat = corr(Activitymatrix');
% correlationmat(logical(eye(size(correlationmat))))= 0;

% num = 3:3:15;
% list1 = sprintf('N1/P%d ',num);
% list2 = sprintf('N2/P%d ',num);
% list = [list1,list2];
% labels = regexp(list,' ','split');
% labels = labels(1:10);
correlationmat(logical(eye(size(correlationmat))))= 0;
figure;
% Plot the pairwise correlations
subplot(2,3,1);
colormap(cmap2);
imagesc(correlationmat);
freezeColors;
if exist('position','var') == 0 || labelopt == 1
    set(gca, 'FontSize', 12)
    xlabel('Neuron #','FontSize', 14);
    ylabel('Neuron #','FontSize', 14);
else
    set(gca,'xtick', tick, 'xticklabel', labels, 'FontSize', 6);
%     set(gca,'xtick', 3:3:30, 'xticklabel', labels) %temp
    ax = gca;
    ax.XTickLabelRotation = 90;
    set(gca,'ytick', tick, 'yticklabel', labels);
%     set(gca,'ytick', 3:3:30, 'yticklabel', labels) %temp
    xlabel('Position (um)');
    ylabel('Position (um)');
%     xlabel('Neuron/Presentation')%temp
%     ylabel('Neuron/Presentation')%temp
end

% line([0 112],[28.5 28.5],'Color','r','LineStyle','--','LineWidth',0.05);
% line([0 112],[56.5 56.5],'Color','r','LineStyle','--','LineWidth',0.05); 
% line([0 112],[84.5 84.5],'Color','r','LineStyle','--','LineWidth',0.05); 
% 
% line([28.5 28.5],[0 112],'Color','r','LineStyle','--','LineWidth',0.05);
% line([56.5 56.5],[0 112],'Color','r','LineStyle','--','LineWidth',0.05);
% line([84.5 84.5],[0 112],'Color','r','LineStyle','--','LineWidth',0.05);

tickpref;


% Plot the spike train matrix
subplot(2,3,[2 3]); 
colormap(cmap2)
zSpikeCount = zscore(Activitymatrix');
imagesc(zSpikeCount');
freezeColors;
set(gca, 'FontSize', 12)
% 
% line([startx startx+500],[28.5 28.5],'Color','r','LineStyle','--','LineWidth',0.05);
% line([startx startx+500],[56.5 56.5],'Color','r','LineStyle','--','LineWidth',0.05);
% line([startx startx+500],[84.5 84.5],'Color','r','LineStyle','--','LineWidth',0.05);

% imagesc(Activitymatrix);
xlim([startx startx+500]);
if exist('position','var') == 0 || labelopt == 1
    ylabel('Neuron #', 'FontSize', 14);
else
    set(gca,'ytick', tick, 'yticklabel', labels, 'FontSize', 6);
    ylabel('Position (um)');
%     set(gca,'ytick', 3:3:30, 'yticklabel', labels) %temp
%     ylabel('Neuron/Presentation') %temp
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
xlabel('Time(ms)', 'FontSize', 14);
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;



% Plot the cell assemblies / independent components
subplot(2,3,4);
colormap(cmap1)
imagesc(Patterns);
freezeColors;
xlabel('Ensemble #');
if exist('position','var') == 0 || labelopt == 1
    ylabel('Neuron #');
else
    set(gca,'ytick', tick, 'yticklabel', labels, 'FontSize', 6);
    ylabel('Position (um)');
%     set(gca,'ytick', 3:3:30, 'yticklabel', labels) %temp
%     ylabel('Neuron/Presentation') %temp
    
end
tickpref;
set(gca,'xtick', 1:size(Patterns,2), 'xticklabel', 1:size(Patterns,2));



% Plot the activities of the cell assemblies / time course of cell assembly activity
subplot(2,3,[5 6]);
plot(Activities');
xlim([startx startx+500]);
xlabel('Time (ms)');
tickpref;
[nr, nc] = size(Activities');
if ( nc > 0 )
    for i = 1:nc
        leg{i} = num2str(i);
    end
    legend(leg,'Location','Best');
end
xtick = get(gca,'xtick');
xticklabel = dt * xtick;
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
set(gcf,'position', [432 334 1121 620]);


% Place a superior title over all the plots
exp = exp_site_nedata.exp;
site = exp_site_nedata.site;
stim = exp_site_nedata.stim;

set(0,'defaulttextinterpreter','none')
% suptitle(sprintf('Cell assembly data of sparsely firing neurons, 100 ms time bin\nNeurons 43 (FR ~ 1.2 Hz) and 52 (FR ~ 2.5 Hz), 50 um apart'))
% suptitle(sprintf('Cell assembly data of densely firing neurons, 100 ms time bin\nNeurons 23 (FR ~ 21 Hz) and 53 (FR ~ 7 Hz), 150 um apart'))
suptitle(sprintf('Exp %s site%.0f %s %.1f ms', ...
    exp, site, stim, ...
    dt));

if nargout >= 2

    varargout{1} = correlationmat;
    varargout{2} = zSpikeCount';
    
elseif nargout >= 1
    
    varargout{1} = correlationmat;
end

% for i = 1:size(Patterns,2)
%     figure;
%     stem(Patterns(:,i));
%     title(sprintf('Ind Comp #%.0f', i));
% end




return;



