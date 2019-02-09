function ne_plot_neuron_NE_firing_rates(neuronrate, NErate, alpha)


figure;
hold on

if iscell(neuronrate)
    try
        neuronrate = cell2mat(neuronrate);
    catch
        neuronrate = cell2mat(neuronrate');
    end
end

if iscell(NErate)
    try
        NErate = cell2mat(NErate);
    catch
        NErate = cell2mat(NErate');
    end
end

meanNErate = mean(NErate, 2);
semNErate = std(NErate,[],2)./sqrt(size(NErate,2));

meanneuronrate = mean(neuronrate);
semneuronrate = std(neuronrate)./sqrt(length(neuronrate));

h = fill([alpha fliplr(alpha)],[meanneuronrate-semneuronrate * ones(1,length(alpha))...
    meanneuronrate+semneuronrate * ones(1,length(alpha))],[200/255 200/255 200/255]);
set(h,'EdgeColor','None');

line([min(alpha) max(alpha)], [meanneuronrate meanneuronrate],'Color', 'k');

errorbar(alpha, meanNErate, semNErate);

xlabel('Threshold percentile')
ylabel('Firing rate (Hz)')
xlim([min(alpha)-0.5 max(alpha)+0.5])