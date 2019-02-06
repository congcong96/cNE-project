function [sil, proj, catmat] = ne_plot_waveform_pca_clusters(exp_site_nedata, spk, edges, rawfolder, varargin)

ip = inputParser;
addRequired(ip, 'exp_site_nedata', @isstruct);
addRequired(ip, 'spk', @isstruct);
addRequired(ip, 'edges', @isvector);
addRequired(ip, 'rawfolder', @ischar);
addParameter(ip, 'groupopt', 'i/a', @(x) strcmp(x, 'i/a') || strcmp(x, 'shared'))
addParameter(ip, 'neuopt', [], @isvector);
addParameter(ip, 'plotopt', 0, @(x) x == 0 || x == 1)
parse(ip, exp_site_nedata, spk, edges, rawfolder, varargin{:})

fn = fieldnames(ip.Results);

for i = 1:length(fn) 
    eval(sprintf('%s = ip.Results.%s;', fn{i}, fn{i}));
end

xbounds = -20:40;

switch groupopt
    case 'i/a' %cNE-i vs cNE-a
        
        spksubset = ne_get_ensemble_spktrain_subset(exp_site_nedata);
        sil = cell(length(spksubset), 1);
        catmat = cell(length(spksubset), 1);
        proj = cell(length(spksubset), 1);
        if plotopt
            f = 0;    
            figure('Position', [349 86 1104 898]);
        end

        for i = 1:length(spksubset)

            fprintf('\nProcessing subset %d of %d',i,length(spksubset))
            w = spksubset(i).spktrain_w_NE;
            wo = spksubset(i).spktrain_wo_NE;
            n = spksubset(i).neuron;

            w_spktime = ne_get_spktime_subset_from_spktrain_subset(spk.spk(n).spiketimes,...
                w(1:end-1), edges);
            wo_spktime = ne_get_spktime_subset_from_spktrain_subset(spk.spk(n).spiketimes,...
                wo(1:end-1), edges);

            w_wavemat = plot_neuron_waveform(spk, rawfolder, n, 'plotopt', 'none', 'xbounds', xbounds, 'subsetspk', w_spktime);
            wo_wavemat = plot_neuron_waveform(spk, rawfolder, n, 'plotopt', 'none', 'xbounds', xbounds, 'subsetspk', wo_spktime);

            wavemat = [wo_wavemat; w_wavemat];
            catmat{i} = [zeros(length(wo_spktime), 1); ones(length(w_spktime), 1)];

            autocorr = corr(wavemat);

            [E,~] = eig(autocorr);

            PCs = E(:,1:2);
            proj{i} = wavemat*PCs;  
            sil{i} = silhouette(proj{i}, catmat{i});

            
                        
            if plotopt
                
                f = f+1;
               
                if f == 5
                    f = 1;
                    figure('Position', [349 86 1104 898]);
                end

                subplot(2,2,f)
                hold on
                scatter(proj{i}(1:length(wo_spktime),1), proj{i}(1:length(wo_spktime),2),10,'filled');
                scatter(proj{i}(length(wo_spktime)+1:end,1), proj{i}(length(wo_spktime)+1:end,2),10,'filled');
                xlabel('PC1')
                ylabel('PC2')
                title(sprintf('Spike waveforms in PC space for cNE-i and cNE-a spikes for neuron %d',...
                    n))
                hold off
                legend('cNE-i','cNE-a')
                tickpref;
                
            end



        end
    
        fprintf('\n');
    
    case 'shared'
        data = ne_plot_shared_neurons_sta(exp_site_nedata, 'stasimopt', 1, 'neuopt', neuopt);
        if isempty(data)
            sil = [];
            proj = [];
            catmat = [];
            return
        end
        spkcell = {data.NEexclusive};
        spkcount = {data.NEexclusive_count};
        n = [data.neuron];
        if plotopt
            f = 0;
            figure('Position', [349 86 1104 898]);
        end

        
        sil = cell(length(spkcell), 1);
        catmat = cell(length(spkcell), 1);
        proj = cell(length(spkcell), 1);

        for i = 1:length(spkcell)
            fprintf('\nProcessing neuron %d of %d',i,length(spkcell))        
            s = size(spkcell{i},1);
            
            spktime = cell(s, 1);
            idx = 1;
            wavemat = zeros(sum(spkcount{i}), length(xbounds));

            
            for ii = 1:s

                spktime{ii} = ne_get_spktime_subset_from_spktrain_subset(spk.spk(n(i)).spiketimes,...
                    spkcell{i}(ii,1:end-1), edges);
                
                numspikes = spkcount{i}(ii);

                wavemat(idx:idx+numspikes-1,:) = plot_neuron_waveform(spk, rawfolder, n(i), 'plotopt', 'none', 'xbounds', xbounds, 'subsetspk', spktime{ii});
                catmat{i}(idx:idx+numspikes-1) = ii*ones(length(spktime{ii}),1);
                
                idx = idx+numspikes;

            end


            autocorr = corr(wavemat);

            [E,~] = eig(autocorr);

            PCs = E(:,1:2);
            proj{i} = wavemat*PCs;
            sil{i} = silhouette(proj{i}, catmat{i});

            
            if plotopt
                
                f = f+1;
               
                if f == 5
                    f = 1;
                    figure('Position', [349 86 1104 898]);
                end

                subplot(2,2,f)
                hold on
                endpt = 0;
                startpt = 1;
                for iii = 1:s
                    endpt = endpt + spkcount{i}(iii); 
                    scatter(proj{i}(startpt:endpt,1),proj{i}(startpt:endpt,2),10,'filled')
                    startpt = startpt + spkcount{i}(iii);
                end
                xlabel('PC1')
                ylabel('PC2')
                title(sprintf('Spikes in PC space with different NEs for neuron %d',...
                    n(i)))
                hold off
                leg = arrayfun(@(x) sprintf('cNE #%d', x), data(i).NEs, 'UniformOutput', 0); 
                legend(leg)
                tickpref;
                
            end


        end

        fprintf('\n');
    
end


