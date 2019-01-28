function [actcorr, trainingact, testact] = ne_predict_subset_of_NE_activity_with_ICweights(exp_site_nedata, varargin)

p = inputParser;
addRequired(p, 'exp_site_nedata', @isstruct)
% addParameter(p, 'dft', 20, @isscalar)
% addParameter(p, 'subset', 0.25, @(x) x>=0 && x<=1)
% addParameter(p, 'sortopt', 0, @(x) x == 0 || x == 1)
addParameter(p, 'plotopt', 0, @(x) x == 0 || x == 1)
addParameter(p, 'plottime', 0, @(x) x >= 0)
parse(p, exp_site_nedata, varargin{:});

exp_site_nedata = p.Results.exp_site_nedata;
% dft = p.Results.dft;
% subset = p.Results.subset;
% sortopt = p.Results.sortopt;
plotopt = p.Results.plotopt;
plottime = p.Results.plottime;

dt = 20;
% niter = 100;

spktrain = exp_site_nedata.nedata.spktrain;
trainlen = floor(size(spktrain,2) / 2);


trainingtrain = spktrain(:,1:trainlen);
testtrain = spktrain(:,trainlen+ 1 : end);


testpat = assembly_patterns(testtrain);
trainingpat = assembly_patterns(trainingtrain);

close all

% if sortopt == 1
% patcorr = abs(corr(trainingpat, testpat));
% [maxrow, colidx] = max(patcorr, [], 2);
% maxidx = [(1:length(maxrow))' colidx];
% % threshpatcorr = patcorr(:, sigcols);
% [~, maxidx] = max(threshpatcorr, [], 1);
% sortedtestpat = testpat(:,sigcols);
% sortedtrainingpat = trainingpat(:,maxidx);
% else %sortopt == 0
% sortedtestpat = testpat;
% sortedtrainingpat = trainingpat;
% end

testact = assembly_activity(testpat, testtrain);
trainingact = assembly_activity(trainingpat, testtrain);

% set baseline activity to 0
% realsubtract = mode(realact, 2);
% testsubtract = mode(testact, 2);
% 
% for i = 1:size(realact, 1)
%     realact(i,:) = realact(i,:) - realsubtract(i);
%     testact(i,:) = testact(i,:) - testsubtract(i);
% end

trainlen = size(testtrain,2);
actcorr = corr(trainingact', testact');

if plotopt == 1
    
%     figure;
%     c = 1;

    [maxrow, colidx] = max(actcorr, [], 2);
    maxidx = [(1:length(maxrow))' colidx];

    for i = 1:size(trainingact, 1)
        figure;
%         subplot(2,2,c)
        hold on

        plot((1:trainlen),trainingact(maxidx(i,1),:))
        plot((5:trainlen+4),testact(maxidx(i,2),:))
        
        xlim([plottime plottime + 1000])        
        xtick = get(gca,'xtick');
        xticklabel = dt * 0.5 / 1000 * xtick;
        set(gca, 'xtick', xtick, 'xticklabel', xticklabel);
        
        
        xlabel('time (s)')
        ylabel('NE activity')
        
%         x = xlim;
        y = ylim;
        
        text(plottime + 500, y(2) - (y(2)-y(1))/2, sprintf('correlation = %.2f', actcorr(maxidx(i,1),maxidx(i,2))))
       
        legend('training', 'test')
        hold off
        tickpref;
        print_mfilename(mfilename)
%         c = c+1;
%         if c == 5
%             figure;
%             c = 1;
%         end
    end
    
end




