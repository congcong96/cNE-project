function varargout = ne_plot_actual_coincidence_vs_model_coincidence(ratio, mtype, statopt)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% Updated 7/21/16 by JS
% Updated 9/7/16 by JS, added choice between mean and median.

if nargin == 1
    mtype = 'median';
    statopt = 0;
elseif nargin == 2
    statopt = 0;
end

% initialize variables
h1 = figure;
h2 = figure;
labels = fieldnames(ratio);
% xlabels = regexprep(labels, '(?<=(\w+))_(?=(ratio))',' ');

% find real ratio
% rridx = find(~cellfun('isempty', strfind(labels, 'real')));
rridx = find(contains(labels, 'real'));
realratio = ratio.(labels{rridx});

% otherratio = rmfield(ratio, labels{rridx});
% otherlabels = labels(~rridx);

c = 1;

for i = 1:length(realratio)
    
    % for each set of N, organize each combination in each row and each
    % manipulation in each column
    mat = [];
    for j = 1:length(labels)
        if size(ratio.(labels{j}){1},2) > 1
            temp1 = cell2mat(ratio.(labels{j})(i,:));
            mat = [mat mean(temp1,2)];
        else
            mat = [mat ratio.(labels{j}){i}(:)];
        end
    end
    mat(isnan(mat)) = 0;
    mat(isinf(mat)) = 0;
    
    % get normalized matrix for coincidence explained
    normmat = bsxfun(@rdivide, mat, mat(:,rridx));
    normmat(isnan(normmat)) = 0;
    normmat(isinf(normmat)) = 0;
    % remove real coincidences with zero values
    zeroidx = normmat(:,rridx) == 0;
    normmat(zeroidx,:) = [];

    
    %% Plot raw data
    set(0, 'currentfigure', h1)
  
    s1 = subplot(2,2,c);
    ax1 = gca;
    hold on

    for k = 1:length(realratio{i})
        plot(1:length(labels) ,mat(k,:),'-o','Color',[0.8,0.8,0.8],'MarkerSize',4)
    end
    
    miu = mean(mat);
    sem = std(mat)./sqrt(size(mat,1));
    errorbar(1:length(labels), miu, sem, 'k', 'LineWidth', 2);
    
    xlim(ax1,[0 length(labels)+ 1])
    tick = 1:length(labels);
    set(ax1,'xtick', tick, 'xticklabel',labels)
    ylabel(ax1,'Coincidence ratio')
    title(s1, sprintf('Coincident activity between %d random sets of %d neurons',...
        length(realratio{i}), i+1))
%     set(gcf, 'Position', [200 100 1500 900])
    
    if length(labels) > 4
        ax1.XTickLabelRotation = 45;
    end
    tickpref;
    hold off    
    
    
    %% Plot normalized data
    set(0, 'currentfigure', h2)
    s2 = figure;
    ax2 = gca;
    hold on
    
    for k = 1:size(normmat,1)
        plot(1:length(labels) ,normmat(k,:),'-o','Color',[0.8,0.8,0.8],'MarkerSize',4)
    end
    
    switch mtype
        case 'median'
            med = median(normmat);
            plot(1:length(labels), med, 'b','LineWidth',1);
            boxinput = normmat(:);
            bsvec = bsxfun(@times, 1:length(labels), ones(size(normmat,1),length(labels)));
            bsvec = bsvec(:);

            boxplot(boxinput, bsvec, 'notch','on','labels',labels,...
                'colors','b', 'symbol','b','outliersize',3, 'widths', 0.5, 'MedianStyle', 'target');
            
        case 'mean'
            miu = mean(normmat);
            sem = std(normmat)./sqrt(size(normmat,1));
            errorbar(1:length(labels), miu, sem, 'k', 'LineWidth', 2);
    end            
    set(ax2, 'tickdir', 'out', 'ticklength', [0.01 0.01]);
    xlim(ax2,[0 length(labels) + 1])
%     ylim(ax2,[0 2])
    tick = 1:length(labels);
    set(ax2,'xtick', tick, 'xticklabel',labels)
    ylabel(ax2,'Normalized coincidence ratio')
    title(sprintf('Coincident activity between %d random sets of %d neurons',...
        size(normmat,1), i+1))
    ylim([0 1.2])
    
    print_mfilename(mfilename)
    
    if length(labels) > 4
        ax2.XTickLabelRotation = 45;
    end
    

    hold off
    
    
    if mod(c, 4) == 0
        set(h1, 'Position', [200 0 1500 900])
        set(h2, 'Position', [200 0 1500 900])
        h1 = figure;
        h2 = figure;
        c = 0;
    end
    c = c+1;
    
    %% Statistics
    if statopt == 1
        
        switch mtype
            case 'median'

                if length(labels) == 2

                    p(i) = signrank(mat(:,1),mat(:,2));
                    str = sprintf('P value for signed-rank test:\n%0.4e',p(i));
                    x = xlim;
                    x = x(2)-((x(2)-x(1))/2.5);
                    y = ylim;
                    y = y(2)-((y(2)-y(1))/5);
                    text(x,y, str)

                else 
                    [~,~,stats] = friedman(mat,1, 'off');
                    figure;
                    hold on
                    p{i} = multcompare(stats,'CType','bonferroni');
                    title(sprintf('Set of %d',i+1))
                    hold off
                end
                
            case 'mean'
                
                
                for j = 1:length(labels)
                    if j ~= rridx
                        [~, temp] = ttest(normmat(:,j), normmat(:,rridx));
                        p{i}(j) = temp * (length(labels) - 1);
                    end
                    
                end
                
        end
                    
                    
        
    end
            
end

if statopt == 1
    varargout{1} = p;
end
    
set(h1, 'Position', [200 0 1500 900])
set(h2, 'Position', [200 0 1500 900])


    

