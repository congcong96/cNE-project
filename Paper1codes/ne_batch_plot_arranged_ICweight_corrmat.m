function [corrmat, varargout] = ne_batch_plot_arranged_ICweight_corrmat(sponfiles, evokedfiles, squareopt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('squareopt','var')
    squareopt = 'none';
else
    color = fifteen_color_blind_palette(10);
end

cmap = brewermap(1000,'reds');
corrmat = cell(length(sponfiles), 1);
maxidx = cell(length(sponfiles), 1);

for i = 1:length(sponfiles)
    
    load(sponfiles{i})
    s_pat = exp_site_nedata.nedata.Patterns;
    load(evokedfiles{i})
    e_pat = exp_site_nedata.nedata.Patterns;
%     e_pat = e_pat(:,1:5); %temp 1/18/19
      
    figure; hold on
    colormap(cmap);
    corrmat{i} = abs(corr(s_pat, e_pat));    
    imagesc(corrmat{i});
      
    xlabel('Evoked cNEs')
    ylabel('Spontaneous cNEs')
    c = colorbar;
    c.TickDirection = 'out';
    ylabel(c, 'Correlation value')
    tickpref;
    print_mfilename(mfilename)
    
    switch squareopt
        case 'max_column'
            
            [~, maxidx{i}] = max(corrmat{i}, [], 1);
            xidx = 0.5:size(corrmat{i},2) - 0.5;
            yidx = maxidx{i} - 0.5;
            
            for j = 1:length(xidx)
                rectangle('Position', [xidx(j) yidx(j) 1 1], 'EdgeColor', color, 'LineWidth', 2)
            end
            
        case 'max_row'
            
            [~, maxidx{i}] = max(corrmat{i}, [], 2);
            xidx = maxidx{i} - 0.5;
            yidx = 0.5:size(corrmat{i}, 2) - 0.5;
            
            for j = 1:length(xidx)
                rectangle('Position', [xidx(j) yidx(j) 1 1], 'EdgeColor', color, 'LineWidth', 2)
            end
                                    
        case 'none'
            continue
                     
                     
    end
    
    xlim([0.49 size(corrmat{i},2)+0.51])
    ylim([0.49 size(corrmat{i},1)+0.51])
    box on

end

varargout{1} = maxidx;
