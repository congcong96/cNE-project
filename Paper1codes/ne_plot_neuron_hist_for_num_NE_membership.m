function [countvecratio, p] = ne_plot_neuron_hist_for_num_NE_membership(files, maxcount)

if ~exist('maxcount','var')
    maxcount = 4;
end

countvec_sig = zeros(1,maxcount + 2);
countvec_nonsig = zeros(1,maxcount + 2);
% count = [];
for i = 1:length(files)
    
    load(files{i}, 'exp_site_nedata')
    members = sort(cell2mat(exp_site_nedata.nedata.NEmembers));
    sig_neurons = find(exp_site_nedata.nedata.sig_neuron_sta);
    nonsig_neurons = find(~exp_site_nedata.nedata.sig_neuron_sta);
    % get neurons part of at least one NE
    [len, val] = rude(members);
    % get neurons not part of any NE
%     uniquemem = unique(members);
    
    [~,sig_idx] = intersect(val, sig_neurons);
    countvec_sig(1) = countvec_sig(1) + length(setdiff(sig_neurons, val));
    [N,~] = histcounts(len(sig_idx),0.5:maxcount+0.5);
    countvec_sig(2:end-1) = countvec_sig(2:end-1) + N;
    countvec_sig(end) = countvec_sig(end) + sum(len(sig_idx) > maxcount);
    
    [~,nonsig_idx] = intersect(val, nonsig_neurons);
    countvec_nonsig(1) = countvec_nonsig(1) + length(setdiff(nonsig_neurons, val));
    [N,~] = histcounts(len(nonsig_idx),0.5:maxcount+0.5);
    countvec_nonsig(2:end-1) = countvec_nonsig(2:end-1) + N;
    countvec_nonsig(end) = countvec_nonsig(end) + sum(len(nonsig_idx) > maxcount);
    
% 
%     countvec_all(1) = size(exp_site_nedata.nedata.spktrain,1) - length(uniquemem) + countvec_all(1);    
%     [N,~] = histcounts(len,0.5:maxcount+0.5);
%     countvec_all(2:end-1) = countvec_all(2:end-1) + N;
%     countvec_all(end) = countvec_all(end) + sum(len > maxcount);

end

bins = 0:maxcount+1;
countvec_all = countvec_sig + countvec_nonsig;

figure;
bar(bins, countvec_all)
% histogram(categorical(count))

xlabel('number of NEs each neuron is a member of')
ylabel('count')

totalcount = sum(countvec_all);

countvecratio = countvec_all ./ sum(countvec_all);
tickpref;
print_mfilename(mfilename)

figure;
bar(bins, countvecratio)
% histogram(categorical(count),'Normalization','probability')

xlabel('number of NEs each neuron is a member of')
ylabel('ratio')

x = 2.8;
y = 0.6;

text(x, y, sprintf('n = %d', totalcount))

tickpref;
print_mfilename(mfilename)



color = eight_color_blind_palette('vermillion', 'bluegreen');

% figure;
% b = bar(bins, [countvec_all' countvec_sig' countvec_nonsig']);
% b(1).FaceColor = color(1,:);
% b(2).FaceColor = color(2,:);
% b(3).FaceColor = color(3,:);
% legend(sprintf('all, n = %d', sum(countvec_all)), sprintf('sig, n = %d',...
%     sum(countvec_sig)), sprintf('nonsig, n = %d', sum(countvec_nonsig)))
% xlabel('cNE membership')
% ylabel('Count')
% tickpref;
% print_mfilename(mfilename)
% 

 

figure;
b = bar(bins, [countvec_sig'./sum(countvec_sig) countvec_nonsig'./sum(countvec_nonsig)]);
for i = 1:size(color, 1)
    b(i).FaceColor = color(i,:);
end
legend(sprintf('sig, n = %d',...
    sum(countvec_sig)), sprintf('nonsig, n = %d', sum(countvec_nonsig)))
xlabel('cNE membership')
ylabel('Ratio')
tickpref;
print_mfilename(mfilename)

exp_proportion = countvec_nonsig ./ sum(countvec_nonsig);
sig_expected = exp_proportion .* sum(countvec_sig);
% [~, p_sig, st_sig] = chi2gof(bins, 'Ctrs', bins, 'Frequency', countvec_sig, 'Expected', sig_expected);


exp_proportion = countvec_sig ./ sum(countvec_sig);
nonsig_expected = exp_proportion .* sum(countvec_nonsig);
% [~, p_nonsig, st_nonsig] = chi2gof(bins, 'Ctrs', bins, 'Frequency', countvec_sig, 'Expected', nonsig_expected);

sample_sig = cell2mat(arrayfun(@(x,y) repmat(x, [1 y]), bins, countvec_sig, 'UniformOutput', 0));
sample_nonsig = cell2mat(arrayfun(@(x,y) repmat(x, [1 y]), bins, countvec_nonsig, 'UniformOutput', 0));

[~, p] = kstest2(sample_sig, sample_nonsig);


% color = eight_color_blind_palette('vermillion', 'bluegreen');
% 
% cat = [repmat({sprintf('sig, n = %d',length(allneurons))},...
%     length(allneurons), 1); repmat({sprintf('non-sig, n = %d',...
%     length(NEneurons))}, length(NEneurons), 1)];
% g1 = gramm('x', [allneurons; NEneurons], 'color', cat);
% g1.stat_bin('edges', edges, 'normalization', 'probability', 'geom', 'line');
% g1.set_names('x', 'Depth (\mum)', 'y', 'Ratio', 'color', 'Depth');
% g1.set_color_options('n_color', 2, 'n_lightness', 1, 'map', color);
% g1.set_order_options('x', 0);
% g1.axe_property('TickDir','out', 'ticklength', [0.01 0.01]);
% figure;
% g1.draw();
% print_mfilename(mfilename); 
