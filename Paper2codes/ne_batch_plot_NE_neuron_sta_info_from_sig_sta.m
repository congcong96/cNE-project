function p = ne_batch_plot_NE_neuron_sta_info_from_sig_sta(sig_sta, plottype, sigopt, paramopt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('paramopt', 'var')
    paramopt = 1;
end
% 
% if ~exist('logopt', 'var')
%     logopt = 1;
% end

prevfilename = [];

NEinfo = cell(length(sig_sta), 1);
neuroninfo = cell(length(sig_sta), 1);
NEptd = cell(length(sig_sta), 1);
neuronptd = cell(length(sig_sta), 1);
NEmi = cell(length(sig_sta), 1);
neuronmi = cell(length(sig_sta), 1);
NEri = cell(length(sig_sta), 1);
neuronri = cell(length(sig_sta), 1);

for i = 1:length(sig_sta)
    
    filename = sig_sta(i).filename;
    
    if ~strcmp(prevfilename, filename)
        load(filename, 'NEneuroninfoall', 'NEneuronstats')
    end
    
    prevfilename = filename;
    
%     load(files{i}, 'NEneuroninfo')
%     load(files{i}, 'NEneuroninfoall')
%     NEneuroninfo = NEneuroninfoall;
    
    if isempty(NEneuroninfoall)
        continue
    else
        
%         NEneuroninfoall(cellfun('isempty', {NEneuroninfoall.min_spikes})) = [];
%         NEneuronstats(cellfun('isempty', {NEneuronstats.min_spikes})) = [];

        
        % find relevant entries in NEneuroninfo struct array
        NE = sig_sta(i).NE;
        
        NEidx = NE == [NEneuroninfoall.NE];
        tempstructMI = NEneuroninfoall(NEidx);
        tempstructstats = NEneuronstats(NEidx);
        
        if sigopt == 1
            sigidx = sig_sta(i).sig_neurons;
        elseif sigopt == 0
            sigidx = ~sig_sta(i).sig_neurons;
        elseif sigopt == -1
            sigidx = true(length(tempstructstats), 1);
        end
        
        if isempty(tempstructMI)
            continue
        end
        
        % remove instances where spk/event count does not cross threshold

        switch plottype
            
            case 'all' % plot cNE vs all neurons
                
                NEinfo{i} = [tempstructMI(sigidx).NE_info_extrap];
                neuroninfo{i} = [tempstructMI(sigidx).neuron_info_extrap];
                
                NEptd{i} = cellfun(@nanmean, {tempstructstats(sigidx).ptd_NE});
                neuronptd{i} = cellfun(@nanmean, {tempstructstats(sigidx).ptd_neuron});
                
                NEmi{i} = cellfun(@nanmean, {tempstructstats(sigidx).moransI_NE});
                neuronmi{i} = cellfun(@nanmean, {tempstructstats(sigidx).moransI_neuron});
                
                NEri{i} = cellfun(@(x) nanmean(x(:)), {tempstructstats(sigidx).reliability_idx_NE});
                neuronri{i} = cellfun(@(x) nanmean(x(:)), {tempstructstats(sigidx).reliability_idx_neuron});
                
%                 if length(NEinfo{i}) ~= length(NEptd{i})
%                     keyboard
%                 end
                
            case 'max' % plot cNE vs best neuron only
                
                diffinfo = [tempstructMI.neuron_info_extrap] - [tempstructMI.NE_info_extrap];
                [~, maxidx] = max(diffinfo);

                NEinfo{i} = [tempstructMI(maxidx).NE_info_extrap];
                neuroninfo{i} = [tempstructMI(maxidx).neuron_info_extrap];
                
                diffptd = cellfun(@mean, {tempstructstats.ptd_neuron}) - cellfun(@mean, {tempstructstats.ptd_NE});
                [~, maxidx] = max(diffptd);
                
                NEptd{i} = mean(tempstructstats(maxidx).ptd_NE);
                neuronptd{i} = mean(tempstructstats(maxidx).ptd_neuron);
                
            case 'mean' % plot cNE vs mean neuron info
                
                NEinfo{i} = nanmean([tempstructMI(sigidx).NE_info_extrap]);
                neuroninfo{i} = nanmean([tempstructMI(sigidx).neuron_info_extrap]);
                
                NEptd{i} = nanmean(cellfun(@nanmean, {tempstructstats(sigidx).ptd_NE}));
                neuronptd{i} = nanmean(cellfun(@nanmean, {tempstructstats(sigidx).ptd_neuron}));
                
                NEmi{i} = nanmean(cellfun(@nanmean, {tempstructstats(sigidx).moransI_NE}));
                neuronmi{i} = nanmean(cellfun(@nanmean, {tempstructstats(sigidx).moransI_neuron}));
                
                NEri{i} = nanmean(cellfun(@(x) nanmean(x(:)), {tempstructstats(sigidx).reliability_idx_NE}));
                neuronri{i} = nanmean(cellfun(@(x) nanmean(x(:)), {tempstructstats(sigidx).reliability_idx_neuron}));
                
            otherwise
                error('Plottype can only be ''all'', ''max'' or ''mean''')
        end
                
    end

end

NEinfo = cell2mat(NEinfo');
neuroninfo = cell2mat(neuroninfo');

NEptd = cell2mat(NEptd');
neuronptd = cell2mat(neuronptd');
nanidx = isnan(NEptd);
NEptd(nanidx) = [];
neuronptd(nanidx) = [];

NEmi = cell2mat(NEmi');
neuronmi = cell2mat(neuronmi');
NEmi(nanidx) = [];
neuronmi(nanidx) = [];

NEri = cell2mat(NEri');
neuronri = cell2mat(neuronri');
NEri(nanidx) = [];
neuronri(nanidx) = [];


figure('Position', [262 84 1330 861]); 
subplot(221); hold on
% if logopt
%     scatter(log10(neuroninfo), log10(NEinfo), 10)
% else
%     scatter(neuroninfo, NEinfo, 10)
% end
scatter(neuroninfo, NEinfo, 6, 'filled')

% if logopt
set(gca, 'xscale','log','yscale','log')
% end
    
h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
tickpref;

% if logopt
%     xlabel('Log(Neuron MI)')
%     ylabel('Log(cNE MI)')
% else
%     xlabel('Neuron MI (bits/spike)')
%     ylabel('cNE MI (bits/event)')
% end
xlabel('Neuron MI (bits/spike)')
ylabel('cNE MI (bits/event)')

if ~paramopt
    p(1) = signrank(NEinfo, neuroninfo);
else
    [~,p(1)] = ttest(NEinfo, neuroninfo);
end

print_n_and_p(3/4, 1/4, length(NEinfo), p(1))
print_mfilename(mfilename);


subplot(222); hold on

scatter(neuronptd, NEptd, 6, 'filled')

h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
tickpref;

xlabel('Neuron PTD')
ylabel('cNE PTD')

if ~paramopt
    p(2) = signrank(NEptd, neuronptd);
else
    [~,p(2)] = ttest(NEptd, neuronptd);
end

print_n_and_p(3/4, 1/4, length(NEptd), p(2))
print_mfilename(mfilename);


subplot(223); hold on
scatter(neuronmi, NEmi, 6, 'filled')

h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
tickpref;

xlabel('Neuron Moran''s I')
ylabel('cNE Moran''s I')

if ~paramopt
    p(3) = signrank(NEmi, neuronmi);
else
    [~,p(3)] = ttest(NEmi, neuronmi);
end

print_n_and_p(3/4, 1/4, length(NEmi), p(3))
print_mfilename(mfilename);

subplot(224); hold on
scatter(neuronri, NEri, 6, 'filled')


h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
tickpref;

xlabel('Neuron reliability index')
ylabel('cNE reliability index')

if ~paramopt
    p(4) = signrank(NEri, neuronri);
else
    [~,p(4)] = ttest(NEri, neuronri);
end

print_n_and_p(3/4, 1/4, length(NEptd), p(4))
print_mfilename(mfilename);

% figure;
% inforatio = NEinfo./neuroninfo;
% histogram(inforatio, 'Normalization', 'probability')
% medval = median(inforatio);
% medad = mad(inforatio, 1);
% y = ylim;
% line([medval medval], [y(1), y(2)], 'LineStyle', '--', 'Color', 'k')
% xlabel('cNE ri / neuron ri')
% ylabel('Ratio')
% tickpref;
% print_mfilename(mfilename);