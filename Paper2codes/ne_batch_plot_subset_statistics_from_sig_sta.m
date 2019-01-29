function ne_batch_plot_subset_statistics_from_sig_sta(sig_sta, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ip = inputParser;
addRequired(ip, 'sig_sta', @isstruct);
addParameter(ip, 'neuronsig', -1, @(x) x == 0 || x == 1 || x == -1);
addParameter(ip, 'NEsig', -1, @(x) x == 0 || x == 1 || x == -1);
parse(ip, sig_sta, varargin{:});

sig_sta = ip.Results.sig_sta;
neuronsig = ip.Results.neuronsig;
NEsig = ip.Results.NEsig;


fn_suffix = {'ptd','moransI','reliability_idx','info_extrap'};
fn_prefix = {'all','w_NE','wo_NE'};

statstruct(length(sig_sta)).temp = [];

prevfile = '';

for i = 1:length(sig_sta)
    
    if ~strcmp(sig_sta(i).filename, prevfile)
        load(sig_sta(i).filename, 'exp_site_nedata', 'subsetstats')
    end
    
    prevfile = sig_sta(i).filename;
    
    neurons = sig_sta(i).neurons;
    NE = sig_sta(i).NE;
    
    subsetneurons = [subsetstats.neuron];
    subsetNEs = [subsetstats.NE];
    ssidx = subsetNEs == NE & logical(sum(cell2mat(arrayfun(@(x) x == subsetneurons, neurons, 'UniformOutput', 0)), 1));
    
    sig_neuron = exp_site_nedata.nedata.sig_neuron_sta;
    sig_NE = exp_site_nedata.nedata.sig_NE_sta;
    
    if neuronsig == -1
        sig_neuron = true(length(sig_neuron), 1);
        neuronidx = 1;
    else
        neuronidx = neuronsig;
    end
    if NEsig == -1
        sig_NE = true(length(sig_NE), 1);
        NEidx = 1;
    else
        NEidx = NEsig;
    end
    
    try
        neuron_num = [subsetstats.neuron]';
        NE_num = [subsetstats.NE]';          
        idx = ssidx' & sig_neuron(neuron_num) == neuronidx & sig_NE(NE_num) == NEidx;
    catch % for penetrations where none of the neurons pass spikecount threshold
        continue
    end
    
    if sum(idx) > 0
        
        for j = 1:length(fn_suffix)
            
            for k = 1:length(fn_prefix)
                
                statstruct(i).([fn_prefix{k} '_' fn_suffix{j}]) = [subsetstats(idx).([fn_prefix{k} '_' fn_suffix{j}])];
                
            end
        end
    end

end

statstruct = rmfield(statstruct, 'temp');
fn = fieldnames(statstruct);
grps = length(fn_prefix);

for i = 1:length(fn)/3    
    
    tempall = cell2mat({statstruct.(fn{(i-1)*grps+1})});
    tempwith = cell2mat({statstruct.(fn{(i-1)*grps+2})});
    tempwithout = cell2mat({statstruct.(fn{(i-1)*grps+3})});

    figure('Position',[324 85 1247 897]);
    
    subplot(2,2,1)
    scatter(tempall, tempwith, 6, 'filled');    
    if contains(fn{(i-1)*grps+1}, 'info') %|| contains(fn{(i-1)*grps+1}, 'ptd')
        set(gca, 'xscale','log','yscale','log')
    end    
    h = refline(1,0);
    h.Color = 'r';
    h.LineStyle = '--';
    tickpref;
    xlabel(sprintf('%s all spikes', fn_suffix{i}));
    ylabel(sprintf('%s cNE-associated spikes', fn_suffix{i}));
    pval = signrank(tempall, tempwith);
    print_n_and_p(1/2, 1/8, length(tempall), pval)

    subplot(2,2,2)
    scatter(tempwithout, tempwith, 6, 'filled');
    if contains(fn{(i-1)*grps+1}, 'info') %|| contains(fn{(i-1)*grps+1}, 'ptd')
        set(gca, 'xscale','log','yscale','log')
    end 
    h = refline(1,0);
    h.Color = 'r';
    h.LineStyle = '--';
    tickpref;
    xlabel(sprintf('%s cNE-independent spikes', fn_suffix{i}));
    ylabel(sprintf('%s cNE-associated spikes', fn_suffix{i}));
    pval = signrank(tempwithout, tempwith);
    print_n_and_p(1/2, 1/8, length(tempwithout), pval)

    subplot(2,2,3)
    scatter(tempwithout, tempall, 6, 'filled');
    if contains(fn{(i-1)*grps+1}, 'info') %|| contains(fn{(i-1)*grps+1}, 'ptd')
        set(gca, 'xscale','log','yscale','log')
    end 
    h = refline(1,0);
    h.Color = 'r';
    h.LineStyle = '--';
    tickpref;
    xlabel(sprintf('%s cNE-independent spikes', fn_suffix{i}));
    ylabel(sprintf('%s all spikes', fn_suffix{i}));
    pval = signrank(tempwithout, tempall);
    print_n_and_p(1/2, 1/8, length(tempwithout), pval)


end

