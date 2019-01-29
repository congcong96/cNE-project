function p = ne_batch_plot_sig_vs_nonsig_NE_depth_span_size(nefiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sigdepth = cell(length(nefiles), 1);
nonsigdepth = cell(length(nefiles), 1);
sigspan = cell(length(nefiles), 1);
nonsigspan = cell(length(nefiles), 1);
sigsize = cell(length(nefiles), 1);
nonsigsize = cell(length(nefiles), 1);
sigproportion = cell(length(nefiles), 1);
nonsigproportion = cell(length(nefiles), 1);
sigperccoin_NE = cell(length(nefiles), 1);
nonsigperccoin_NE = cell(length(nefiles), 1);
sigperccoin_neuron = cell(length(nefiles), 1);
nonsigperccoin_neuron = cell(length(nefiles), 1);
sigNEsigneuron_perccoin = cell(length(nefiles), 1);
nonsigNEnonsigneuron_perccoin = cell(length(nefiles), 1);
sigNEnonsigneuron_perccoin = cell(length(nefiles), 1);
nonsigNEsigneuron_perccoin = cell(length(nefiles), 1);
sigpwd = cell(length(nefiles), 1);
nonsigpwd = cell(length(nefiles), 1);

color = eight_color_blind_palette('orange', 'blue');

for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    
    depth = cellfun(@(x) x(2), exp_site_nedata.nedata.position);
    
    sigdepth{i} = cellfun(@(x) median(depth(x)), exp_site_nedata.nedata.NEmembers...
        (exp_site_nedata.nedata.sig_NE_sta));
    
    nonsigdepth{i} = cellfun(@(x) median(depth(x)), exp_site_nedata.nedata.NEmembers...
        (~exp_site_nedata.nedata.sig_NE_sta));
    
    sigspan{i} = cellfun(@(x) max(depth(x)) - min(depth(x)), exp_site_nedata.nedata.NEmembers...
        (exp_site_nedata.nedata.sig_NE_sta));
    
    nonsigspan{i} = cellfun(@(x) max(depth(x)) - min(depth(x)), exp_site_nedata.nedata.NEmembers...
        (~exp_site_nedata.nedata.sig_NE_sta));
    
    sigsize{i} = cellfun('length', exp_site_nedata.nedata.NEmembers(...
        exp_site_nedata.nedata.sig_NE_sta));
    
    nonsigsize{i} = cellfun('length', exp_site_nedata.nedata.NEmembers(...
        ~exp_site_nedata.nedata.sig_NE_sta));
    
    sigproportion{i} = sigsize{i} ./ length(exp_site_nedata.nedata.position);
    
    nonsigproportion{i} = nonsigsize{i} ./ length(exp_site_nedata.nedata.position);
    
    NEmemtrain = ne_get_neuronal_ensemble_spktrain(exp_site_nedata,...
        'threshalpha', 99.5, 'method', 'repeat', 'memneuopt', 'raw');
    
    perc_coin = cellfun(@(x) sum(x(2:end, logical(x(1,:))), 2)...
        ./ sum(x(2:end, :), 2), NEmemtrain, 'UniformOutput',0);
    
    sigperccoin_NE{i} = cell2mat(perc_coin(exp_site_nedata.nedata.sig_NE_sta));
    
    nonsigperccoin_NE{i} = cell2mat(perc_coin(~exp_site_nedata.nedata.sig_NE_sta));
    
    sigidx_NE = cell2mat(arrayfun(@(x, y) repmat(x, [y 1]), exp_site_nedata.nedata.sig_NE_sta,...
        cellfun('length', exp_site_nedata.nedata.NEmembers), 'UniformOutput', 0));    
    perc_coin = cell2mat(perc_coin);
    sigidx_neuron = exp_site_nedata.nedata.sig_neuron_sta(cell2mat(exp_site_nedata.nedata.NEmembers));
    
    sigperccoin_neuron{i} = perc_coin(sigidx_neuron);
    
    nonsigperccoin_neuron{i} = perc_coin(~sigidx_neuron);
    
    sigNEsigneuron_perccoin{i} = perc_coin(sigidx_neuron & sigidx_NE);
    sigNEnonsigneuron_perccoin{i} = perc_coin(~sigidx_neuron & sigidx_NE);
    nonsigNEnonsigneuron_perccoin{i} = perc_coin(~sigidx_neuron & ~sigidx_NE);
    nonsigNEsigneuron_perccoin{i} = perc_coin(sigidx_neuron & ~sigidx_NE);
    
    neurons = exp_site_nedata.nedata.NEmembers;

    tempcell = cell(length(neurons), 1);
    for k = 1:length(neurons)
        if length(neurons{k}) <= 1
            continue
        end
        comb = nchoosek(neurons{k}, 2);
        tempcell{k} = ne_calc_interneuronal_distances(comb,...
            cell2mat(exp_site_nedata.nedata.position'), 'distance');
    end
    sigpwd{i} = cell2mat(tempcell(exp_site_nedata.nedata.sig_NE_sta));
    nonsigpwd{i} = cell2mat(tempcell(~exp_site_nedata.nedata.sig_NE_sta));


end


depthstruct.sig = cell2mat(sigdepth);
depthstruct.nonsig = cell2mat(nonsigdepth);
figure('Position', [680 478 372 500]);
p{1} = plot_plotspread_and_boxplot(depthstruct, 'ranksum', color, -1);
ylabel('cNE depth (um)')
axis ij
print_mfilename(mfilename);

spanstruct.sig = cell2mat(sigspan);
spanstruct.nonsig = cell2mat(nonsigspan);
figure('Position', [680 478 372 500]);
p{2} = plot_plotspread_and_boxplot(spanstruct, 'ranksum', color, 1);
ylabel('cNE span (um)')
print_mfilename(mfilename);

sizestruct.sig = cell2mat(sigsize);
sizestruct.nonsig = cell2mat(nonsigsize);
figure('Position', [680 478 372 500]);
p{3} = plot_plotspread_and_boxplot(sizestruct, 'ranksum', color, 1);
ylabel('cNE size')
print_mfilename(mfilename);

propstruct.sig = cell2mat(sigproportion);
propstruct.nonsig = cell2mat(nonsigproportion);
figure('Position', [680 478 372 500]);
p{4} = plot_plotspread_and_boxplot(propstruct, 'ranksum', color, 1);
ylabel('cNE size (ratio)')
print_mfilename(mfilename);


perccoinstructNE.sig_NE = cell2mat(sigperccoin_NE);
perccoinstructNE.nonsig_NE = cell2mat(nonsigperccoin_NE);
figure('Position', [680 478 372 500]);
p{5} = plot_plotspread_and_boxplot(perccoinstructNE, 'ranksum', color, 1);
ylabel('Percentage spike coincidence')
print_mfilename(mfilename);

perccoinstructneuron.sig_neuron = cell2mat(sigperccoin_neuron);
perccoinstructneuron.nonsig_neuron = cell2mat(nonsigperccoin_neuron);
figure('Position', [680 478 372 500]);
p{6} = plot_plotspread_and_boxplot(perccoinstructneuron, 'ranksum', color, 1);
ylabel('Percentage spike coincidence')
print_mfilename(mfilename);

pwdstruct.sig = cell2mat(sigpwd);
pwdstruct.nonsig = cell2mat(nonsigpwd);
figure('Position', [680 478 372 500]);
p{7} = plot_plotspread_and_boxplot(pwdstruct, 'ranksum', color, 1);
ylabel('Pairwise distance (um)')
print_mfilename(mfilename);


color = eight_color_blind_palette('vermillion','redpurple','bluegreen','skyblue');
perccoinstruct.facilitative = cell2mat(sigNEsigneuron_perccoin);
perccoinstruct.constructive = cell2mat(sigNEnonsigneuron_perccoin);
perccoinstruct.independent = cell2mat(nonsigNEnonsigneuron_perccoin);
perccoinstruct.destructive = cell2mat(nonsigNEsigneuron_perccoin);
figure('Position', [680 478 372*2 500]);
p{8} = plot_plotspread_and_boxplot(perccoinstruct, 'kruskalwallis', color, 1);
ylabel('Percentage spike coincidence')
print_mfilename(mfilename);
