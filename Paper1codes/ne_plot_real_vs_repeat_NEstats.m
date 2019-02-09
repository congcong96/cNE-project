function p = ne_plot_real_vs_repeat_NEstats(NEstats)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Plot numNE
realnumNE = NEstats.realnumNE;
repnumNE = NEstats.repnumNE;

figure;
hold on

plotSpread(realnumNE, 'distributionColors', [.8 .8 .8], 'xValues', 1, 'spreadWidth', 0.5);
plotSpread(repnumNE, 'distributionColors', [.8 .8 .8], 'xValues', 2, 'spreadWidth', 0.5);

mnrealnumNE = mean(realnumNE);
semrealnumNE = std(realnumNE);%./sqrt(length(realnumNE));

mnrepnumNE = mean(repnumNE);
semrepnumNE = std(realnumNE);%./sqrt(length(repnumNE));

errorbar(1:2, [mnrealnumNE mnrepnumNE], [semrealnumNE semrepnumNE], 'kx')
set(gca,'xtick', 1:2, 'xticklabel',{'real', 'repeat shuffled'})

ymax = ceil(max([realnumNE;repnumNE]))+1;
ymin = floor(min([realnumNE;repnumNE]))-1;

ylim([ymin ymax])

ylabel('Number of NEs')
tickpref;

p.numNE = ranksum(realnumNE, repnumNE);

hold off

%Plot NEsize
realNEsize = NEstats.realNEsize;
repNEsize = NEstats.repNEsize;

figure;
hold on

plotSpread(realNEsize, 'distributionColors', [.8 .8 .8], 'xValues', 1, 'spreadWidth', 0.5);
plotSpread(repNEsize, 'distributionColors', [.8 .8 .8], 'xValues', 2, 'spreadWidth', 0.5);

mnrealNEsize = mean(realNEsize);
semrealNEsize = std(realNEsize);%./sqrt(length(realNEsize));

mnrepNEsize = mean(repNEsize);
semrepNEsize = std(repNEsize);%./sqrt(length(repNEsize));

errorbar(1:2, [mnrealNEsize mnrepNEsize], [semrealNEsize semrepNEsize], 'kx')
set(gca,'xtick', 1:2, 'xticklabel',{'real', 'repeat shuffled'})

ymax = round(max([realNEsize;repNEsize]))+1;
ymin = round(min([realNEsize;repNEsize]))-1;

ylim([ymin ymax])

ylabel('Mean NE size')
tickpref;

p.NEsize = ranksum(realNEsize, repNEsize);

hold off


end

