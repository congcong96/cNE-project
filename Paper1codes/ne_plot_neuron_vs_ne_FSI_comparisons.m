function [FSI_paired, p] = ne_plot_neuron_vs_ne_FSI_comparisons(FSI_ne, FSI_neuron, FSIthresh, comparisontype, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch comparisontype
    
    case 'unpaired'

        if iscell(FSI_ne)
            FSI_ne = cell2mat(FSI_ne);
        end

        if iscell(FSI_neuron)
            FSI_neuron = cell2mat(FSI_neuron);
        end

        for i = 1:length(FSIthresh)

            figure; hold on
            NE = FSI_ne(:,i);
            neuron = FSI_neuron(:,i);

            plotSpread(NE, 'distributionColors', [.8 .8 .8], 'xValues',...
                0, 'spreadWidth', 1)

            plotSpread(neuron, 'distributionColors', [.8 .8 .8], 'xValues',...
                2, 'spreadWidth', 1)

            boxinput = [NE;neuron];
            bsvec = [zeros(length(NE),1);ones(length(neuron),1)];
            labels = {'cNEs','neurons'};

            boxplot(boxinput, bsvec, 'notch','on','labels', labels, 'colors','b',...
                'symbol','b','outliersize',3, 'widths', 0.5, 'position', [0 2])

            ylabel('Ratio of distribution')
            title(sprintf('CDF value at %.1f', FSIthresh(i)))

            [~, p(1)] = ttest2(NE, neuron);
            p(2) = ranksum(NE,neuron);

            text(1,mean(NE),sprintf('p-param = %.3f\np-nonparam = %.3f',p(1), p(2)));
        %     xlim([-1 1])
        end
        
    case 'paired'
       
        files = varargin{1};
        
        c = 1;
        
        for k = 1:length(files)
            
            load(files{k}, 'exp_site_nedata')
            NEmembers = exp_site_nedata.nedata.NEmembers;
            numNEmembers = sum(cellfun('length', NEmembers));
%             FSI_paired = zeros(numNEmembers, 2);          

            for i = 1:length(NEmembers)

                for j = 1:length(NEmembers{i})
                    
                    for ii = 1:length(FSIthresh)
                        
                        FSI_paired{ii}(c,1) = FSI_ne{k}(i,ii);
                        FSI_paired{ii}(c,2) = FSI_neuron{k}(NEmembers{i}(j), ii);
                        
                    end
                c = c+1;
                end
            end
            
        end
        
        
        for i = 1:length(FSI_paired)
            
            figure; hold on
            scatter(FSI_paired{i}(:,1), FSI_paired{i}(:,2));
            h = refline(1,0);
            h.Color = 'r';
            h.LineStyle = '--';
            title(sprintf('FSI = %.1f', FSIthresh(i)));
            
            [~,p(1,i)] = ttest(FSI_paired{i}(:,1), FSI_paired{i}(:,2));
            p(2,i) = signrank(FSI_paired{i}(:,1), FSI_paired{i}(:,2));
            
        end
        
        
end

