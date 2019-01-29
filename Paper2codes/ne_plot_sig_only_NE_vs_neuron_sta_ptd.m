function [p, stastruct] = ne_plot_sig_only_NE_vs_neuron_sta_ptd(NEneuronsta, logopt)

% Plot peak-trough difference for significant cNE-neuron pairs only.

% NEneuronsta = cell(length(files), 1);
% 
% for i = 1:length(files)
%     
%     fprintf('\nProcessing %s...\n', files{i})
%     
%     load(files{i}, 'exp_site_nedata')
%     
%     NEneuronsta{i} = ne_calc_NE_vs_neuron_ptd(exp_site_nedata);
% 
% end

NEneuronsta = NEneuronsta(~cellfun('isempty', NEneuronsta));

stastruct = cell2mat(NEneuronsta');
minspikes = [stastruct.min_spikes];

ptdNE = cellfun(@mean, {stastruct.ptd_NE}) ./minspikes;
ptdneuron = cellfun(@mean, {stastruct.ptd_neuron}) ./minspikes;


figure; hold on
if logopt
    scatter(log10(ptdneuron), log10(ptdNE), 10)
    xlabel('Log(Neuron PTD)')
    ylabel('Log(cNE PTD)')
else
    scatter(ptdneuron, ptdNE, 10)
    xlabel('Neuron PTD')
    ylabel('cNE PTD')
end
h = refline(1,0);
h.Color = [.7 .7 .7];
h.LineStyle = '--';
tickpref;



p(1) = signrank(ptdNE, ptdneuron);
[~,p(2)] = ttest(ptdNE, ptdneuron);

text(2*10^4,10^4,sprintf('n = %d\np < 0.001', length(ptdNE)));

print_mfilename(mfilename);

figure;
ptdratio = ptdNE./ptdneuron;
histogram(ptdratio, 'Normalization', 'probability')
medval = median(ptdratio);
medad = mad(ptdratio, 1);
y = ylim;
line([medval medval], [y(1), y(2)], 'LineStyle', '--', 'Color', 'k')
xlabel('cNE PTD / neuron PTD')
ylabel('Ratio')
tickpref;
print_mfilename(mfilename);



