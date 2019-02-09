function ne_plot_pfr_comparison_examples(spktrain, spktrainPred, binsize, numsamp, settime)
%UNTITLED10 Summary of this function goes here
%   binsize is in seconds
if nargin == 3
    numsamp = 4;
    settime = [];
elseif nargin == 4
    settime = [];
end

fr = mean(spktrain);
frPred = mean(spktrainPred);

if isempty(settime)    
    mindiff = 0;
    while mindiff < 600
        randidx = randsample(1:length(fr)-600, numsamp);
        mindiff = min(diff(sort(randidx)));
    end
    randidx = sort(randidx);    
else
    assert(length(settime) == numsamp)
    settime = settime/binsize;
    randidx = settime;
end
    
figure;
for i = 1:numsamp
    subplot(2,2,i)
    hold on
    plot(randidx(i)*binsize:binsize:(randidx(i)+500)*binsize, fr(randidx(i):randidx(i) + 500))
    plot(randidx(i)*binsize:binsize:(randidx(i)+500)*binsize, frPred(randidx(i):randidx(i) + 500))
    hold off
    xlim([randidx(i)*binsize (randidx(i)+500)*binsize])
    xtick = get(gca,'xtick');
    set(gca, 'xtick', xtick, 'xticklabel', 0:5)
    ylabel('Population firing rate')
    xlabel('Time / s')
    if i == numsamp
        legend('real spiketrain', 'surrogate spiketrain')
    end
    tickpref;
end

set(gcf, 'Position', [200 200 1500 750])

