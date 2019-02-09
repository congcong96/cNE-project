function varargout = ne_plot_exp_site_nedata_fig3(exp_site_nedata, startx)
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
end


nedata = exp_site_nedata.nedata;
Activitymatrix = nedata.spktrain;


Patterns = nedata.Patterns;
Activities = nedata.Activities;
NEmembers = nedata.NEmembers;


try
    dt = nedata.df/2; % spike train bin size, in ms
catch
    dt = nedata.binsize; 
end


% Plot the pairwise correlations
correlationmat = corr(Activitymatrix');
correlationmat(logical(eye(size(correlationmat))))= 0;

cmapcorr = disproportionate_divergent_colormap(min(correlationmat(:)),...
    max(correlationmat(:)), 1000, 'rdbu', 'brange', [0 0.5]);


figure;
colormap(cmapcorr);
imagesc(correlationmat);
colorbar;
set(gca, 'FontSize', 12)
xlabel('Neuron #','FontSize', 14);
ylabel('Neuron #','FontSize', 14);

tickpref;


% Plot the spike train matrix
figure; 

[~, ~] = plotSpikeRaster(logical(Activitymatrix), 'PlotType','vertline');
% cmaprast = greyscale_zscored_raster(Activitymatrix);
% colormap(cmaprast)
% zSpikeCount = zscore(Activitymatrix');
% imagesc(zSpikeCount');

xlim([startx startx+200]);

set(gca, 'FontSize', 12)
ylabel('Neuron #', 'FontSize', 14);

xtick = get(gca,'xtick');
xticklabel = 0:0.5:2;
xlabel('Time(s)', 'FontSize', 14);
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
tickpref;
set(gcf,'position', [432 334 1121 620]);
box on


% Plot the cell assemblies / independent components

% Flip ICs if they are dominantly negative
markers = cell(length(NEmembers),1);

for i = 1:size(Patterns,2)
    patsum = sum(Patterns(:,i));
    markers{i} = [ones(length(NEmembers{i}),1)*i NEmembers{i}];
    
    if patsum < 0
        Patterns(:,i) = -Patterns(:,i);
    end
end

figure;
markersmat = cell2mat(markers);
cmapic = disproportionate_divergent_colormap(min(Patterns(:)),...
    max(Patterns(:)), 1000, 'rdbu');

colormap(cmapic);
hold on;
imagesc(Patterns);
s = scatter(markersmat(:,1), markersmat(:,2), 20, [0 .8 0],'filled');
s.Marker = 'square';
box on
axis ij
colorbar;
xlabel('Ensemble #');
ylabel('Neuron #');

xlim([0.5 length(NEmembers)+0.5])
ylim([0.5 size(Patterns,1)+0.5])

tickpref;
set(gca,'xtick', 1:size(Patterns,2), 'xticklabel', 1:size(Patterns,2));



% Plot the activities of the cell assemblies / time course of cell assembly activity
figure;
hold on
linecolors = distinguishable_colors(size(Activities,1));
leg = cell(size(Activities,1),1);
for i = 1:size(Activities,1)
    plot(Activities(i,:), 'Color', linecolors(i,:));
    leg{i} = num2str(i);
end
    
xlim([startx startx+200]);
xlabel('Time (s)', 'FontSize', 14);
tickpref;    
legend(leg,'Location','Best');
box on
set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
set(gcf,'position', [432 334 1121 620]);
hold off

Activities = Activities([2 5],:);
n = size(Activities,1);
figure;
for i= 1:n
    subplot(n, 1, i)
    plot(Activities(i,:));
    xlim([startx startx+200]);
    set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
    tickpref;
end

set(gcf,'position', [432 334 1121 620]);


if nargout >= 2

    varargout{1} = correlationmat;
    varargout{2} = zSpikeCount';
    
elseif nargout >= 1
    
    varargout{1} = correlationmat;
end

return;



