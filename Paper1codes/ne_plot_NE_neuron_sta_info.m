function p = ne_plot_NE_neuron_sta_info(files, plottype, logopt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

NEinfo = [];
neuroninfo = [];

for i = 1:length(files)
    
    load(files{i}, 'NEneuroninfo')
%     load(files{i}, 'NEneuroninfoall')
%     NEneuroninfo = NEneuroninfoall;
    
    if isempty(NEneuroninfo)
        continue
    else
        
        % remove instances where spk/event count does not cross threshold
        NEneuroninfo(cellfun('isempty', {NEneuroninfo.min_spikes})) = [];

        switch plottype
            
            case 'all' % plot cNE vs all neurons
                NEinfo = [NEinfo NEneuroninfo.NE_info_extrap];
                neuroninfo = [neuroninfo NEneuroninfo.neuron_info_extrap];
                
            case 'max' % plot cNE vs best neuron only
                NEs = [NEneuroninfo.NE];
                uniqueNEs = unique(NEs);
                diffinfo = [NEneuroninfo.neuron_info_extrap] - [NEneuroninfo.NE_info_extrap];

                idx = zeros(length(uniqueNEs), 1);

                for j = 1:length(uniqueNEs)
                    tempidx = find(NEs == uniqueNEs(j));
                    [~, maxidx] = max(diffinfo(tempidx));
                    idx(j) = tempidx(maxidx);
                end
                
                tempNEinfo = [NEneuroninfo.NE_info_extrap];
                tempneuroninfo = [NEneuroninfo.neuron_info_extrap];

                NEinfo = [NEinfo tempNEinfo(idx)];
                neuroninfo = [neuroninfo tempneuroninfo(idx)];
                
            case 'mean' % plot cNE vs mean neuron info
                
                NEs = [NEneuroninfo.NE];
                uniqueNEs = unique(NEs);
                allNEinfo = [NEneuroninfo.NE_info_extrap];
                allneuroninfo = [NEneuroninfo.neuron_info_extrap];

                tempNEinfo = zeros(1,length(uniqueNEs));
                tempneuroninfo = zeros(1, length(uniqueNEs));
                for j = 1:length(uniqueNEs)
                    tempidx = NEs == uniqueNEs(j);
                    tempNEinfo(j) = mean(allNEinfo(tempidx));
                    tempneuroninfo(j) = mean(allneuroninfo(tempidx));
                end

                NEinfo = [NEinfo tempNEinfo];
                neuroninfo = [neuroninfo tempneuroninfo];
            
            otherwise
                error('Plottype can only be ''all'', ''max'' or ''mean''')
        end
                
    end

end


figure; hold on
if logopt
    scatter(log10(neuroninfo), log10(NEinfo), 10)
else
    scatter(neuroninfo, NEinfo, 10)
end
h = refline(1,0);
h.Color = [.7 .7 .7];
h.LineStyle = '--';
tickpref;

if logopt
    xlabel('Log(Neuron MI)')
    ylabel('Log(cNE MI)')
else
    xlabel('Neuron MI (bits/spike)')
    ylabel('cNE MI (bits/event)')
end

p(1) = signrank(NEinfo, neuroninfo);
[~,p(2)] = ttest(NEinfo, neuroninfo);
x = xlim;
y = ylim;

if p < 0.001
    text(3/4*(x(2)-x(1))+x(1),y(1)+(y(2)-y(1))/4,sprintf('n = %d\np < 0.001', length(NEinfo)));
else
    text(3/4*(x(2)-x(1))+x(1),y(1)+(y(2)-y(1))/4,sprintf('n = %d\np = %.2f', length(NEinfo), p(1)));
end

print_mfilename(mfilename);


figure;
inforatio = NEinfo./neuroninfo;
histogram(inforatio, 'Normalization', 'probability')
medval = median(inforatio);
medad = mad(inforatio, 1);
y = ylim;
line([medval medval], [y(1), y(2)], 'LineStyle', '--', 'Color', 'k')
xlabel('cNE MI / neuron MI')
ylabel('Ratio')
tickpref;
print_mfilename(mfilename);