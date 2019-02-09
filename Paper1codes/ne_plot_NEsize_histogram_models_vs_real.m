function ne_plot_NEsize_histogram_models_vs_real(NEsize, models)

if ischar(models) && ~strcmp(models, 'all')
    models = {models};
end

fn = fieldnames(NEsize);
realidx = strcmp(fn, 'real');
realstats = NEsize.(fn{realidx});
modelstats = rmfield(NEsize, fn{realidx});

if ischar(models) && strcmp(models,'all')
    
    models = fieldnames(modelstats);

    figure;
    hold on
    box on
    for i = 1:length(models)
        histogram(modelstats.(models{i}))
    end
    
elseif iscell(models)
    figure;
    hold on
    box on
    for i = 1:length(models)
        temp = modelstats.(models{i});
        maxmax = max(temp);
%         if strcmp(models{i}, 'prePFRRF')
%             histogram(modelstats.(models{i}), 'Normalization', 'probability')
%         else
            edges = 0:1:maxmax;
            histogram(modelstats.(models{i}), edges, 'Normalization','probability')
%         end
    end
end

y = ylim;

line([realstats realstats],[y(1) y(2)], 'Color','r','LineStyle','--')

% xlim([-0.5 realstats + 0.5])

text(realstats/2, y(2)/2, 'n = 200');

legend([models; 'real data'], 'Location', 'best')
xlabel('mean NE size')
ylabel('ratio')
tickpref;
hold off
    
    