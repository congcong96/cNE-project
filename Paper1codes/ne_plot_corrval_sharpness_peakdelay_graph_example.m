function ne_plot_corrval_sharpness_peakdelay_graph_example(exp_site_nedata, allpwc, neuron)

% Plot graphical representation of correlation values between spiketrains,
% pwc sharpness, and pwc peak delay


% get correlation values
corrmat = corr(exp_site_nedata.nedata.spktrain');
corrvals = corrmat(:, neuron);
corrvals(neuron) = [];

% get pwc involving one neuron only
pwcpairs = {allpwc.pairs}';
pwcidx = cellfun(@(x) sum(x == neuron) ~= 0, pwcpairs); 
pwc = allpwc(pwcidx);

% get pwc sharpness
[~, ~, interpvals, ~, pwcidx2] = calc_pwc_sharpness_cdf(pwc, 40, 0.5, 0);
% interpvals(isnan(interpvals)) = 0;

% get pwc peakdelay
pd = [pwc.peakdelay]';
peakdelay = abs(pd(pwcidx2));

% get relevant correlation values
corrs = corrvals(pwcidx2);

% get neuron numbers
neuronnum = 1:size(exp_site_nedata.nedata.spktrain,1);
neuronnum(neuronnum == neuron) = [];


A = zeros(size(corrmat));
Aidx = sort([ones(length(interpvals),1)*neuron neuronnum(pwcidx2)'],2);
Aidxcell = num2cell(Aidx,1);
Aidx = sub2ind(size(corrmat), Aidxcell{:});
A(Aidx) = interpvals;

% % define edges
% t = neuronnum(pwcidx2);
% s = ones(length(t),1) * neuron;
% 
% % fill in rest of neurons if not present
% if max(t) ~= size(corrmat,1) && s(1) ~= size(corrmat,1)
%     t = [t size(corrmat,1)];
%     s = [s; neuron];
%     interpvals = [interpvals; 0];
% end
    
% get graph
G = graph(A, 'upper');

LWidths = 10*G.Edges.Weight/max(G.Edges.Weight);
LWidths(isnan(LWidths)) = 0.001;

% plot graph with edge width representing sharpness of pwc and edge label
% representing pwc peak delay
figure;
h = plot(G, 'LineWidth', LWidths, 'EdgeLabel', peakdelay);

% get coordinates of equal points on circumference of circle with radius 1
angles = 0 : 2*pi/length(corrvals) : 2*pi;
angles = angles(1:end-1);
xvals = cos(angles);
yvals = sin(angles);

% reorganize points for reference neuron to be in the middle and the rest
% of the neurons to be in a circle around the ref neuron
h.XData = [xvals(1:neuron-1) 0 xvals(neuron:end)];
h.YData = [yvals(1:neuron-1) 0 yvals(neuron:end)];

% change node color to black for clarity
h.NodeColor = 'k';

% get correlation colors
[cmapcorr,~,~] = brewermap(1000,'reds');
cmapcorr(sum(cmapcorr,2) == 3, :) = [];

% assign each corr value to color
edges = linspace(min(corrs), max(corrs), size(cmapcorr,1) + 1);
[~,~,bin] = histcounts(corrs, edges);
edgecolor = cmapcorr(bin,:);

% change colors of edges according to corr values
for i = 1:length(pwcidx2)
    highlight(h, neuron, neuronnum(pwcidx2(i)), 'EdgeColor', edgecolor(i,:))
end

set(gcf, 'Position', [200 100 1200 800])

figure;
colormap(cmapcorr)
imagesc(corrs)
colorbar;


end

