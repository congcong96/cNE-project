function ne_plot_CI_vs_mean_vs_num_neurons(files)

% sanity check function to test if CI for determination of which neurons
% belong in which NEs makes sense

% Written 10/26/16 by JS

L = length(files);

numneurons = zeros(L,1);
CI = zeros(L,1);
ICweights = cell(L,1);
ICweightsmn = zeros(L,1);
ICweightsmed = zeros(L,1);

for i = 1:length(files)
    
    load(files{i})
    numneurons(i) = size(exp_site_nedata.nedata.spktrain,1);
    CI(i) = mean(abs(exp_site_nedata.nedata.CI));
    ICweights{i} = abs(exp_site_nedata.nedata.Patterns(:));
    ICweightsmn(i) = mean(ICweights{i});
    ICweightsmed(i) = median(ICweights{i});
end


figure;
subplot(2,2,1)
hold on
scatter(ICweightsmn, ICweightsmed)
h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
xlabel('abs IC weight mean')
ylabel('abs IC weight median')

subplot(2,2,2)
scatter(numneurons, ICweightsmed);
xlabel('number of neurons')
ylabel('abs IC weight median')

subplot(2,2,3)
scatter(numneurons, CI);
xlabel('number of neurons')
ylabel('abs mean threshold')

subplot(2,2,4)
scatter(CI, ICweightsmed);
xlabel('abs mean threshold')
ylabel('abs IC weight median')

return