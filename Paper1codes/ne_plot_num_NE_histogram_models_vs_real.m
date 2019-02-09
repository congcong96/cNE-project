function ne_plot_num_NE_histogram_models_vs_real(NEstats, models)

if ischar(models) && ~strcmp(models, 'all')
    models = {models};
end

fn = fieldnames(NEstats);
realidx = strcmp(fn, 'real');
realstats = NEstats.(fn{realidx});
modelstats = rmfield(NEstats, fn{realidx});

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
        histogram(modelstats.(models{i}), 'Normalization','probability')
    end
end

y = ylim;

line([realstats realstats],[y(1) y(2)], 'Color','r','LineStyle','--')

xlim([-0.5 realstats + 0.5])

text(realstats/2, y(2)/2, 'n = 200');

legend([models; 'real data'], 'Location', 'best')
xlabel('number of NEs')
ylabel('count')
tickpref;
hold off
    
    