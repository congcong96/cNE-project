function sil = ne_plot_waveform_pca_clusters(exp_site_nedata, spk, edges, type, neuopt)

if ~exist('neuopt', 'var')
    neuopt = [];
end

switch type
    case 'assembly'
        
        spksubset = ne_get_ensemble_spktrain_subset(exp_site_nedata);
        f = 0;

        for i = 1:length(spksubset)

            fprintf('\nProcessing subset %d of %d',i,length(spksubset))
            w = spksubset(i).spktrain_w_NE;
            wo = spksubset(i).spktrain_wo_NE;
            n = spksubset(i).neuron;

            w_spktime = ca_get_spktime_subset_from_spktrain_subset(spk(n).spiketimes,...
                w, exp_site_nedata.nedata.edges);
            wo_spktime = ca_get_spktime_subset_from_spktrain_subset(spk(n).spiketimes,...
                wo, exp_site_nedata.nedata.edges);

            w_wavemat = plot_neuron_waveform(spk, n, 'none', -20:80, w_spktime);
            wo_wavemat = plot_neuron_waveform(spk, n, 'none', -20:80, wo_spktime);

            w_size = size(w_wavemat,1);

            wavemat = [w_wavemat'; wo_wavemat'];

            autocorr = corr(wavemat);

            [E,~] = eig(autocorr);

            PCs = E(:,1:2);
            proj = wavemat*PCs;

            f = f+1;
            if mod(f,4) == 1
                f = 1;
                figure;
            end


            subplot(2,2,f)
            hold on
            scatter(proj(1:w_size,1),proj(1:w_size,2),10,'r','filled')
            scatter(proj(w_size+1:end,1),proj(w_size+1:end,2),10,'b','filled')
            xlabel('PC1')
            ylabel('PC2')
            legend('withCA','withoutCA')
            title(sprintf('Spikes in PC space with CA and without CA for neuron %d assembly %d',...
                n,spksubset(i).ensemble))
            hold off


        end

        fprintf('\n');
    
    case 'shared'
        data = ne_plot_shared_neurons_sta(exp_site_nedata, 'stasimopt', 1, 'neuopt', neuopt);
        if isempty(data)
            sil = [];
            return
        end
        spkcell = {data.NEexclusive};
        spkcount = {data.NEexclusive_count};
        n = [data.neuron];
        f = 0;
        
    sil = zeros(length(spkcell),1);
        
    for i = 1:length(spkcell)
        fprintf('\nProcessing neuron %d of %d',i,length(spkcell))        
        s = size(spkcell{i},1);
        wavemat = [];
        catmat = [];
        
        for ii = 1:s
            
            spktime = ca_get_spktime_subset_from_spktrain_subset(spk.spk(n(i)).spiketimes,...
                spkcell{i}(ii,1:end-1), edges);
            
            wave = plot_neuron_waveform(spk, n(i), 'none', -20:80, spktime);
            wavemat = [wavemat; wave];
            catmat = [catmat; ii*ones(length(spktime),1)];
            
        end
        
        sil(i) = mean(silhouette(wavemat, catmat));

        autocorr = corr(wavemat);

        [E,~] = eig(autocorr);

        PCs = E(:,1:2);
        proj = wavemat*PCs;

        f = f+1;
        if mod(f,4) == 1
            f = 1;
            figure;
        end


        subplot(2,2,f)
        hold on
        endpt = 0;
        startpt = 1;
        for iii = 1:s
            endpt = endpt + spkcount{i}(iii); 
            scatter(proj(startpt:endpt,1),proj(startpt:endpt,2),10,'filled')
            startpt = startpt + spkcount{i}(iii);
        end
        xlabel('PC1')
        ylabel('PC2')
        title(sprintf('Spikes in PC space with different NEs for neuron %d',...
            n(i)))
        hold off
        legend(sprintf('NE #%d', data(i).NEs(1)), sprintf('NE #%d', data(i).NEs(2)))
        tickpref;


    end

    fprintf('\n');
    
end
