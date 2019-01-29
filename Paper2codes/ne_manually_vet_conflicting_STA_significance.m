function [NEvet, neuronvet] = ne_manually_vet_conflicting_STA_significance(sig_sta1, sig_sta2)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

prompt = 'Non-significant(0), significant(1) or ambiguous(2)?';

%% Plot cNEs

NEidx = find([sig_sta1.sig_NE] ~= [sig_sta2.sig_NE]);

filenames = {sig_sta1(NEidx).filename};

NEvet(length(filenames)).filename = [];
[NEvet.filename] = filenames{:};

figure('Position', [2084 359 976 458]);
cmap = flipud(brewermap(1000, 'rdbu'));
colormap(cmap);

for i = 1:length(NEidx)
    
    NEvet(i).filename = sig_sta1(NEidx(i)).filename;
    NEvet(i).NE = sig_sta1(NEidx(i)).NE;
    NEvet(i).sigNE1 = sig_sta1(NEidx(i)).sig_NE;
    NEvet(i).sigNE2 = sig_sta2(NEidx(i)).sig_NE;
    
    subplot(1,2,1)
    quick_plot_sta(sig_sta1(NEidx(i)).NE_sta);
    
    subplot(1,2,2)
    quick_plot_sta(sig_sta1(NEidx(i)).NE_sig_sta);
    
    suptitle(sprintf('Event count: %d', sig_sta1(NEidx(i)).NE_event_count));
    
    NEvet(i).input_sig = input(prompt);

end

%% Plot neurons

uniqfn = unique({sig_sta1.filename});
c = 1;
neuronvet(length(uniqfn)*10).filename = [];

for i = 1:length(uniqfn)
    
    idx = cellfun(@(x) strcmp(x, uniqfn{i}), {sig_sta1.filename});
    
    [neurons, nidx] = unique(cell2mat({sig_sta1(idx).neurons}'), 'stable');
    
    signeurons1 = cell2mat({sig_sta1(idx).sig_neurons}');
    signeurons1 = signeurons1(nidx);
    signeurons2 = cell2mat({sig_sta2(idx).sig_neurons}');
    signeurons2 = signeurons2(nidx);    
    
    neuron_sta = cell2mat({sig_sta1(idx).neuron_sta}');
    neuron_sta = neuron_sta(nidx, :);
    neuron_sig_sta = cell2mat({sig_sta1(idx).neuron_sig_sta}');
    neuron_sig_sta = neuron_sig_sta(nidx, :);
    
    spikecount = cell2mat({sig_sta1(idx).neuron_spike_count}');
    spikecount = spikecount(nidx, :);
    
    cidx = find(signeurons1 ~= signeurons2);
    
    if isempty(cidx)
        continue
    end
    
    for j = 1:length(cidx)
        
        neuronvet(c).filename = uniqfn{i};
        neuronvet(c).neuron = neurons(cidx(j));
        neuronvet(c).signeuron1 = signeurons1(cidx(j));
        neuronvet(c).signeuron2 = signeurons2(cidx(j));
        
        subplot(1,2,1)
        quick_plot_sta(neuron_sta(cidx(j),:));
        
        subplot(1,2,2)
        quick_plot_sta(neuron_sig_sta(cidx(j),:));
        
        suptitle(sprintf('Spike count: %d', spikecount(cidx(j)))) ;
        
        neuronvet(c).input_sig = input(prompt);
        c = c+1;
        
    end
        
end

neuronvet(cellfun('isempty', {neuronvet.input_sig})) = [];

