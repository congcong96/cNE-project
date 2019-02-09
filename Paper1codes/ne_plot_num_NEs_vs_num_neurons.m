function b = ne_plot_num_NEs_vs_num_neurons(files)

% Plots number of assemblies against number of neurons for a list of files.
% Should only compare across files with similar time bins.
% 
% files: cell array of file names to be calculated.

mat = zeros(length(files),2);
% stim = cell(length(files),1);


for i = 1:length(files)
    
    load(files{i})
    mat(i,:) = size(exp_site_nedata.nedata.Patterns);
%     stim{i} = regexp(files{i}, '(?<=(db-))rn\d{1,3}(?=(-fs))','match','once');
    
end

exp = '(?<=(-ne-))\w+(?=(dft.mat))';
dft = str2double(regexp(files{1},exp,'match','once'));


% color code different stims
% uniqstim = unique(stim);
% markertype = {'s','^','o','d'};
figure;
hold on
% for j = 1:length(uniqstim)
%     idx = strcmp(stim, uniqstim{j});
scatter(mat(:,1), mat(:,2),'filled')%,markertype{j})
% end
% figure;
% hold on
% scatter(mat(:,1),mat(:,2));
X = [ones(size(mat,1),1) mat(:,1)];
Y = mat(:,2);
[b, ~,~,~,stats] = regress(Y,X);
x = min(mat(:,1))-5:max(mat(:,1))+5;
y = b(1) + b(2)*x;
plot(x,y,'k--');


xlabel('Number of neurons')
ylabel('Number of cNEs')
title(sprintf('Number of cNEs against number of neurons for %d-ms timebins',dft*0.5))

yl = ylim;
xl = xlim;
if stats(3) < 0.001
    text(xl(1)+(xl(2)-xl(1))*0.25, yl(1) + (yl(2)-yl(1))*0.75, ...
        sprintf('n = %d\nR^{2} = %.2f\np < 0.001',length(files), stats(1)))
else
    text(xl(1)+(xl(2)-xl(1))*0.25, yl(1) + (yl(2)-yl(1))*0.75,...
        sprintf('n = %d\nR^{2} = %.2f\np = %.3f', length(files), stats(1), stats(3)))
end

% legend(uniqstim, 'Location', 'best')
tickpref;
print_mfilename(mfilename);

end

