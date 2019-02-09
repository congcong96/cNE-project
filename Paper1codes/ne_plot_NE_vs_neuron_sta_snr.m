function p = ne_plot_NE_vs_neuron_sta_snr(files)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

NEneuronsta = cell(length(files), 1);

for i = 1:length(files)
    
    fprintf('\nProcessing %s...\n', files{i})
    
    load(files{i}, 'exp_site_nedata')
    
    NEneuronsta{i} = ne_plot_NE_vs_neuron_sta(exp_site_nedata, 0);
    


end


stastruct = cell2mat(NEneuronsta');
minspikes = [stastruct.min_spikes];

staNE = {stastruct.sta_NE};
staNE = staNE(~cellfun('isempty', staNE));
NEmat = cell2mat(cellfun(@(x) x(:), staNE, 'UniformOutput', 0));
maxNEmat = max(NEmat, [], 1);
minNEmat = min(NEmat, [], 1);
snrNE = (maxNEmat - minNEmat)./minspikes;

staneuron = {stastruct.sta_neuron};
staneuron = staneuron(~cellfun('isempty', staneuron));
neuronmat = cell2mat(cellfun(@(x) x(:), staneuron, 'UniformOutput', 0));
maxneuronmat = max(neuronmat, [], 1);
minneuronmat = min(neuronmat, [], 1);
snrneuron = (maxneuronmat - minneuronmat)./minspikes;


figure; hold on
scatter(snrneuron, snrNE, 10)
h = refline(1,0);
h.Color = [.7 .7 .7];
h.LineStyle = '--';
tickpref;


xlabel('Neuron max - min')
ylabel('cNE max - min')

p(1) = signrank(snrNE, snrneuron);
[~,p(2)] = ttest(snrNE, snrneuron);

text(2*10^4,10^4,sprintf('n = %d\np < 0.001', length(snrNE)));

print_mfilename(mfilename);

figure;
snrratio = snrNE./snrneuron;
histogram(snrratio, 'Normalization', 'probability')
medval = median(snrratio);
medad = mad(snrratio, 1);
y = ylim;
line([medval medval], [y(1), y(2)], 'LineStyle', '--', 'Color', 'k')
xlabel('cNE MI / neuron MI')
ylabel('Ratio')
tickpref;
print_mfilename(mfilename);



