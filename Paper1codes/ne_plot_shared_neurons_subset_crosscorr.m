function ne_plot_shared_neurons_subset_crosscorr(spkcell)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

x = -100:5:100;

for i = 1:length(spkcell)
%     assert(size(spkcell{i},1) == 2)
    crosscorr = xcorr(spkcell{i}(1,:),spkcell{i}(2,:), 20);
    figure;
    bar(x, crosscorr)
    xlim([-100 100])
    xlabel('Delay (ms)')
    ylabel('Count')
    tickpref;
end

