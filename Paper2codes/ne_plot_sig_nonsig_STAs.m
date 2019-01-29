function ne_plot_sig_nonsig_STAs(exp_site_nedata)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here

sig_NE = [];
nonsig_NE = [];
sig_neuron = [];
nonsig_neuron = [];

nf = exp_site_nedata.nedata.nf;
nlags = exp_site_nedata.nedata.nlags;
nedata = exp_site_nedata.nedata;
sig_stats = nedata.sig_classifying_stats;

NEsig = nedata.sig_NE_sta;
sig_NEsta = nedata.sig_NE_stamat;

temp = sig_NEsta(NEsig,:);
sig_NE = [sig_NE; temp];

temp = sig_NEsta(~NEsig,:);
nonsig_NE = [nonsig_NE; temp];

neuronsig = nedata.sig_neuron_sta;
sig_neuronsta = nedata.sig_stamat;

temp = sig_neuronsta(neuronsig,:);
sig_neuron = [sig_neuron; temp];

temp = sig_neuronsta(~neuronsig,:);
nonsig_neuron = [nonsig_neuron; temp];



cmap = flipud(brewermap(1000, 'RdBu'));

figure('units','normalized','outerposition',[0 0 1 1]); 
colormap(cmap);
c = 1;   
label = find(NEsig);
ptd_pval = sig_stats.NE_ptd_pval(NEsig);
moransI_pval = sig_stats.NE_moransI_pval(NEsig);
relidx_pval = sig_stats.NE_relidx_pval(NEsig);

for i = 1:size(sig_NE, 1)

    subplot(4,4,c)
    imagesc(reshape(sig_NE(i,:), nf, nlags))
    title(sprintf('cNE #%d, %.2f, %.2f, %.2f', label(i), ptd_pval(i), moransI_pval(i), relidx_pval(i)))
    cmax = max(sig_NE(i,:));
    cmin = min(sig_NE(i,:));
    climit = max(abs([cmax cmin]));
    set(gca, 'clim',[-climit climit])
    axis xy
    tickpref;

    if c == 16 && i ~= size(sig_NE, 1)
        suptitle('Significant cNE STAs')
        figure('units','normalized','outerposition',[0 0 1 1]); 
        colormap(cmap);
        c = 1;
    else
        c = c+1;
    end
end
suptitle('Significant cNE STAs')

figure('units','normalized','outerposition',[0 0 1 1]); 
colormap(cmap);
c = 1;
label = find(~NEsig);    
ptd_pval = sig_stats.NE_ptd_pval(~NEsig);
moransI_pval = sig_stats.NE_moransI_pval(~NEsig);
relidx_pval = sig_stats.NE_relidx_pval(~NEsig);

for i = 1:size(nonsig_NE, 1)

    subplot(4,4,c)
    imagesc(reshape(nonsig_NE(i,:), nf, nlags))
    title(sprintf('cNE #%d, %.2f, %.2f, %.2f', label(i), ptd_pval(i), moransI_pval(i), relidx_pval(i)))
    cmax = max(nonsig_NE(i,:));
    cmin = min(nonsig_NE(i,:));
    climit = max(abs([cmax cmin]));
    set(gca, 'clim',[-climit climit])
    axis xy
    tickpref;

    if c == 16 && i ~= size(nonsig_NE, 1)
        suptitle('Non-significant cNE STAs')
        figure('units','normalized','outerposition',[0 0 1 1]); 
        colormap(cmap);
        c = 1;
    else
        c = c+1;
    end
end
suptitle('Non-significant cNE STAs')    
    

figure('units','normalized','outerposition',[0 0 1 1]); 
colormap(cmap);
c = 1;
label = find(neuronsig);    
ptd_pval = sig_stats.neuron_ptd_pval(neuronsig);
moransI_pval = sig_stats.neuron_moransI_pval(neuronsig);
relidx_pval = sig_stats.neuron_relidx_pval(neuronsig);
    
for i = 1:size(sig_neuron, 1)

    subplot(4,4,c)
    imagesc(reshape(sig_neuron(i,:), nf, nlags))
    title(sprintf('Neuron #%d, %.2f, %.2f, %.2f', label(i), ptd_pval(i), moransI_pval(i), relidx_pval(i)))
    cmax = max(sig_neuron(i,:));
    cmin = min(sig_neuron(i,:));
    climit = max(abs([cmax cmin]));
    set(gca, 'clim',[-climit climit])
    axis xy
    tickpref;

    if c == 16 && i ~= size(sig_neuron, 1)
        suptitle('Significant neuron STAs')
        figure('units','normalized','outerposition',[0 0 1 1]); 
        colormap(cmap);
        c = 1;
    else
        c = c+1;
    end
end
suptitle('Significant neuron STAs')

figure('units','normalized','outerposition',[0 0 1 1]); 
colormap(cmap);
c = 1;
label = find(~neuronsig); 
ptd_pval = sig_stats.neuron_ptd_pval(~neuronsig);
moransI_pval = sig_stats.neuron_moransI_pval(~neuronsig);
relidx_pval = sig_stats.neuron_relidx_pval(~neuronsig);
    
for i = 1:size(nonsig_neuron, 1)

    subplot(4,4,c)
    imagesc(reshape(nonsig_neuron(i,:), nf, nlags))
    title(sprintf('Neuron #%d, %.2f, %.2f, %.2f', label(i), ptd_pval(i), moransI_pval(i), relidx_pval(i)))
    cmax = max(nonsig_neuron(i,:));
    cmin = min(nonsig_neuron(i,:));
    climit = max(abs([cmax cmin]));
    set(gca, 'clim',[-climit climit])
    axis xy
    tickpref;

    if c == 16 && i ~= size(nonsig_neuron, 1)
        suptitle('Non-significant neuron STAs')
        figure('units','normalized','outerposition',[0 0 1 1]); 
        colormap(cmap);
        c = 1;
    else
        c = c+1;
    end
end
suptitle('Non-significant neuron STAs')    

    
