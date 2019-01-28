function percent = ne_plot_shared_neurons_sta_corr_pval_hist(nefiles, sigopt)

% if iscell(rcorrval)
%     rcorrmat = cell2mat(rcorrval);
%     scorrcell = cat(1,scorrval{:});
% else
%     rcorrmat = rcorrval;
%     scorrcell = scorrval;
% end

if ~exist('sigopt','var')
    sigopt = 1;
end

corrcell = cell(length(nefiles), 1);
for i = 1:length(nefiles)
    load(nefiles{i}, 'shared_STAcorr', 'exp_site_nedata');
    
    if isempty(shared_STAcorr)
        continue
    end
    
    if sigopt
        neuronnum = [shared_STAcorr.neuron];
        signeuron = exp_site_nedata.nedata.sig_neuron_sta;
        idx = signeuron(neuronnum);
        shared_STAcorr = shared_STAcorr(idx);
    end   
    
    corrcell{i} = shared_STAcorr;
end

corrstruct = horzcat(corrcell{:});
rcorrcell = {corrstruct.real_corrvals};
scorrcell = {corrstruct.surr_corrvals};

percent = zeros(length(rcorrcell),1);
medscorrmat = zeros(length(rcorrcell),1);
rcorrval = zeros(length(rcorrcell),1);

for i = 1:length(rcorrcell)
    
    rcorrval(i) = mean(rcorrcell{i});    
    percent(i) = sum(rcorrval(i) >= [scorrcell{i}; rcorrval(i)]) / (length(scorrcell{i} + 1));
    medscorrmat(i) = median(scorrcell{i});
    
end

figure;
edges = 0:0.05:1;
histogram(percent,edges, 'Normalization', 'probability');
% histogram(percent,edges);
tickpref;
text(0.5,0.3, sprintf('n = %d', length(rcorrcell)));
xlabel('Uniqueness index')
ylabel('Ratio')


figure;
% colormap(flipud(brewermap(251,'OrRd')));
colormap(brewermap(251,'RdYlBu'));


scatter(rcorrval, medscorrmat, 6, percent, 'filled');
h = refline(1,0);
h.Color = [.7 .7 .7];
h.LineStyle = '--';
xlabel('Median sim STA correlation (bits/spike)')
ylabel('Real STA correlation (bits/event)')
tickpref;
cb = colorbar;
ylabel(cb, 'Uniqueness index')
set(cb, 'TickDirection', 'out')
print_mfilename(mfilename)
p = signrank(rcorrval, medscorrmat);
print_n_and_p(3/4, 1/4, length(rcorrval), p);

% count = histcounts(percent,edges);
% 
% x = edges(2:end);
% 
% bar(x, count, 1);