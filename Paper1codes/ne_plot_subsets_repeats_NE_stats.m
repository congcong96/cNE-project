function [p1, p2] = ne_plot_subsets_repeats_NE_stats(files)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

numNEmat = zeros(length(files),2);
NEsizemat = numNEmat;

numNE = figure;
NEsize = figure;

for i = 1:length(files)
    
    load(files{i})
    numNEmat(i,1) = mean(repNE.realnumNEs);
    numNEmat(i,2) = mean(repNE.surrnumNEs);
    NEsizemat(i,1) = mean(repNE.realNEsize);
    NEsizemat(i,2) = mean(repNE.surrNEsize);
    
    set(0, 'currentfigure', numNE);
    hold on
    plot([0.5 1], numNEmat(i,:), 'Color', [.8 .8 .8])
    
    set(0, 'currentfigure', NEsize);
    hold on
    plot([0.5 1], NEsizemat(i,:), 'Color', [.8 .8 .8])
end
    

meannumNE = mean(numNEmat);
semnumNE = std(numNEmat)./sqrt(size(numNEmat,1));
meanNEsize = mean(NEsizemat);
semNEsize = std(NEsizemat)./sqrt(size(NEsizemat,1));


set(0, 'currentfigure', numNE);
hold on
errorbar([0.5 1], meannumNE, semnumNE, 'k');
ylabel('mean number of NEs')
set(gca, 'XTick', [0.5 1], 'XTickLabel', {'real NEs','surrogate NEs'})
[~, p1(1)] = ttest(numNEmat(:,1), numNEmat(:,2));
p1(2) = signrank(numNEmat(:,1), numNEmat(:,2));
tickpref;
print_mfilename(mfilename);
 

set(0, 'currentfigure', NEsize);
hold on
errorbar([0.5 1], meanNEsize, semNEsize, 'k');
ylabel('mean NE size')
set(gca, 'XTick', [0.5 1], 'XTickLabel', {'real NEs','surrogate NEs'})
[~, p2(1)] = ttest(NEsizemat(:,1), NEsizemat(:,2));
p2(2) = signrank(NEsizemat(:,1), NEsizemat(:,2));
tickpref;
print_mfilename(mfilename);


end

