function p = ne_plot_model_vs_real_numNE_NEsize(files, model, param, plotopt)

if ischar(model) && strcmp(model, 'all')
    model = {'preISI','preRF','prePFRRF'};
elseif ischar(model)
    model = {model};
end

meanNEs(length(files)).(model{1}) = [];
realNEs = zeros(length(files), 1);

for i = 1:length(files)
    
    load(files{i})
    
    for j = 1:length(model)
        
        switch param
            case 'numNE'
                meanNEs(i).(model{j}) = mean(numNE.(model{j}));
                realNEs(i) = numNE.real;
            case 'NEsize'
                meanNEs(i).(model{j}) = mean(NEsize.(model{j}));
                realNEs(i) = NEsize.real;
            otherwise
                error('Enter either numNE or NEsize for param')
        end
        
    end
end

    
switch plotopt
    case 'separate'

        for k = 1:length(model)

            figure;
            hold on
            
            for ii = 1:length(files)
                plot(1:2, [meanNEs(ii).(model{k}) realNEs(ii)],'-o', 'Color',[.8 .8 .8],'MarkerSize',4)
            end

            mu = mean([[meanNEs.(model{k})]; realNEs], 2);
            sem = std([[meanNEs.(model{k})]; realNEs], 0, 2)./sqrt(length(files));

            errorbar(1:2, mu, sem, 'k', 'LineWidth', 2);

            xlim([0 3])
            
            switch param
                case 'numNE'
                    xlabels = {['mean number of NEs\newline' '(' model{j} ')'],...
                        'number of NEs\newline(real)'};
                    ax = gca;
                    set(ax,'xtick', 1:2, 'xticklabel',xlabels)
                    ylabel('(mean) number of NEs')
                case 'NEsize'
                    xlabels = {['mean NE size\newline' '(' model{j} ')'],...
                        'mean NE size\newline(real)'};
                    ax = gca;
                    set(ax,'xtick', 1:2, 'xticklabel',xlabels)
                    ylabel('mean NE size')
            end
            ax.XTickLabelRotation = 45;
            tickpref;
            p = signrank([meanNEs.(model{k})], realNEs);
            
        end
        print_mfilename(mfilename);


    case 'together'

        figure;
        hold on
        modelNEs = zeros(length(files), length(model));
        
        for ii = 1:length(files)
            modelNEs(ii,:) = cell2mat(struct2cell(meanNEs(ii)));
            plot(1:length(model)+1, [modelNEs(ii,:) realNEs(ii)],'-o', 'Color',[.8 .8 .8],'MarkerSize',4)
        end
        
        mu = mean([modelNEs realNEs]);
        sem = std([modelNEs realNEs])./sqrt(length(files));

        errorbar(1:length(model)+1, mu, sem, 'k', 'LineWidth', 2);

        xlim([0 length(model)+2])
        
        switch param
            case 'numNE'
%                 modellabels = cellfun(@(x) strcat('mean number of NEs\newline',...
%                     '(', x, ')'), model, 'UniformOutput', 0);
                xlabels = [model 'real'];
                ax = gca;
                set(ax,'xtick', 1:length(model)+1, 'xticklabel',xlabels)
                ylabel('(mean) number of NEs')
            case 'NEsize'
%                 modellabels = cellfun(@(x) strcat('mean NE size\newline',...
%                     '(', x, ')'), model, 'UniformOutput', 0);
                xlabels = [model 'real'];
                ax = gca;
                set(ax,'xtick', 1:length(model)+1, 'xticklabel',xlabels)
                ylabel('mean NE size')
        end
%         ax.XTickLabelRotation = 45;
        tickpref;
        
        for k = 1:length(model)
            [~, p(k)] = ttest(modelNEs(:,k), realNEs);
        end
    %Bonferroni correction
    p = p * length(model);
    
    print_mfilename(mfilename);
end


            