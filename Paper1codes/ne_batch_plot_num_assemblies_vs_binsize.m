function assemnum = ne_batch_plot_num_assemblies_vs_binsize(fileid, ratioopt, plotopt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    ratioopt = 1;
    plotopt = 1;
elseif nargin == 2
    plotopt = 1;
end

if length(fileid{1}) > 20
    fileid = cellfun(@(x) x(1:19), fileid, 'UniformOutput',0);
end

assemnum = cell(length(fileid),1);

for i = 1:length(fileid)
    files = gfn([fileid{i} '*dft.mat']);
    stim = unique(regexp(files, '(?<=(db-))\w+(?=(-fs))','match','once'));
    dft = unique(regexp(files, '(?<=(ne-))\d{1,3}(?=(dft.mat))','match','once'));
    assemnum{i} = zeros(length(stim),length(dft));
    
    for j = 1:length(stim)
         [assemnum{i}(j,:), binsize] = ne_get_num_assemblies_vs_binsize(fileid{i}, stim{j}, ratioopt);
    end
end

plotmat = cell2mat(assemnum);

figure;
hold on

if plotopt == 0
    for k = 1:size(plotmat,1)    
        plot(binsize,plotmat(k,:),'o-')    
    end
    
elseif plotopt == 1
    
    for k = 1:size(plotmat,1)
        plot(binsize,plotmat(k,:),'-o','Color',[0.8,0.8,0.8],'MarkerSize',4)
    end
    
    miu = mean(plotmat);
    sem = std(plotmat)./sqrt(size(plotmat,1));
    
%     plot(binsize, miu, '-o', 'Color', [0 0 0], 'MarkerSize',5)
    
    errorbar(binsize, miu, sem, 'k-o')
    title('Number of NEs vs binsize used')
    xlabel('binsize/ms')
    ylabel('normalized number of NEs')
    
end
hold off
tickpref;
print_mfilename(mfilename);