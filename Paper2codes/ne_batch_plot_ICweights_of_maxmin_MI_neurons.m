function p = ne_batch_plot_ICweights_of_maxmin_MI_neurons(nefiles, sigopt, posopt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ICwtMIstruct = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    load(nefiles{i}, 'exp_site_nedata', 'NEneuroninfoall')
    if isempty(NEneuroninfoall)
        continue
    end    
    ICwtMIstruct{i} = ne_calc_ICweights_MIdiff_NEneurons(exp_site_nedata, NEneuroninfoall);
end

ICwtMIstruct = cell2mat(ICwtMIstruct');
ICwtMIstruct(cellfun('isempty', {ICwtMIstruct.NE})) = [];

minMI = zeros(length(ICwtMIstruct),1);
maxMI = zeros(length(ICwtMIstruct),1);

allMIdiff = cell(length(ICwtMIstruct), 1);
allICwt = cell(length(ICwtMIstruct), 1);

for i = 1:length(ICwtMIstruct)
    
    ICwt = ICwtMIstruct(i).ICwt;
    MIdiff = ICwtMIstruct(i).MIdiff;

    if sigopt % takes only significant cNEs         
        if ~ICwtMIstruct(i).sig_NE
            continue
        else
            neuronidx = ICwtMIstruct(i).sig_neurons;
            
            if sum(neuronidx) <= 1 % skip if only 1 significant neuron
                continue
            end
            
            ICwt = ICwt(neuronidx);
            MIdiff = MIdiff(neuronidx);
        end        
    end
    
    if posopt % ignores negative ICwt neurons
        
        posidx = ICwt > 0; 
        if sum(posidx) <= 1 % skip if only 1 neuron with positive IC weight
            continue
        end
        
        ICwt = ICwt(posidx);
        MIdiff = MIdiff(posidx);
    end
            
    [~, minidx] = min(ICwt);
    [~, maxidx] = max(ICwt);
    
    minMI(i) = MIdiff(minidx);
    maxMI(i) = MIdiff(maxidx);
    
    allMIdiff{i} = MIdiff;
    allICwt{i} = ICwt;
    
end

% remove skipped values
zeroidx = minMI == 0 & maxMI == 0;
minMI(zeroidx) = [];
maxMI(zeroidx) = [];


% plot maxMI vs minMI
figure;
scatter(minMI, maxMI, 10, 'filled');
xlabel('MI diff of min IC weight neuron')
ylabel('MI diff of max IC weight neuron')
h1 = refline(1,0);
h1.Color = 'r';
h1.LineStyle = '--';

p(1) = signrank(minMI, maxMI);
[~, p(2)] = ttest(minMI, maxMI);

xlim([-.5, .7])
ylim([-.5, .7])

print_n_and_p(1/5, 4/5, length(minMI), p(1)) 
tickpref;
print_mfilename(mfilename);


% plot MIdiff ICwt scatterplot
figure; hold on
allMIdiff = cell2mat(allMIdiff')';
allICwt = cell2mat(allICwt);

scatter(allICwt, allMIdiff, 10, 'filled');
xlabel('IC weight')
ylabel('MI difference')

[b,~,~,~,stats] = regress(allMIdiff(:),[ones(length(allICwt), 1) allICwt(:)]);

h2 = refline(b(1), b(2));
h2.Color = 'k';
% h2.LineStyle = '--';
print_n_and_p(3/4, 3/4, length(allMIdiff), stats(3)) 
text(0.7, 1.2, sprintf('R^2 = %.2e', stats(1)));
tickpref;
print_mfilename(mfilename);