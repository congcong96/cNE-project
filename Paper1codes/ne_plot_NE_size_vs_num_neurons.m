function [b] = ne_plot_NE_size_vs_num_neurons(files)

% Plots size of NE (mean number of neurons per NE) against number of neurons for a list of files.
% Should only compare across files with similar time bins.
% 
% files: cell array of file names to be calculated.

NEsize = zeros(length(files),1);
numneu = zeros(length(files),1);
% stim = cell(length(files),1);

for i = 1:length(files)
    
    load(files{i})
%     stim{i} = regexp(files{i}, '(?<=(db-))rn\d{1,3}(?=(-fs))','match','once');
    EM = exp_site_nedata.nedata.NEmembers;
    NEsize(i) = mean(cellfun(@length, EM));
    numneu(i) = size(exp_site_nedata.nedata.spktrain,1);
end

% color code different stims
% uniqstim = unique(stim);
% markertype = {'s','^','o','d'};

figure;
hold on
% for j = 1:length(uniqstim)
%     idx = strcmp(stim, uniqstim{j});
scatter(numneu, NEsize, 'filled')%, markertype{j})
% end
% legend(uniqstim, 'Location', 'best')

exp = '(?<=(-ne-))\w+(?=(dft.mat))';
dft = str2double(regexp(files{1},exp,'match','once'));

% figure;
% hold on
% scatter(numneu,NEsize);
X = [ones(length(NEsize),1) numneu];
Y = NEsize;
[b, ~, ~,~, stats] = regress(Y,X);
x = min(numneu)-5:max(numneu)+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');


xlabel('Number of neurons')
ylabel('Mean cNE size')
title(sprintf('Mean cNE size against number of neurons for %d-ms timebins',dft*.5))

yl = ylim;
xl = xlim;

if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.25, yl(1) + (yl(2)-yl(1))*0.75, ...
        sprintf('n = %d\nR^{2} = %.2f\np < 0.001', length(files),stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.25, yl(1) + (yl(2)-yl(1))*0.75,...
        sprintf('n = %d\nR^{2} = %.2f\np = %.3f', length(files), stats(1), stats(3)))
end
tickpref;    
print_mfilename(mfilename);

end

