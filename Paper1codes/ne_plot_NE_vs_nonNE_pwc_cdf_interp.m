function p = ne_plot_NE_vs_nonNE_pwc_cdf_interp(hdNE, hdnNE, indivplot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% if nargin == 2
%     pval = 0.5;
% end

% %eliminate pwcs with meagre counts
% count = {pwcNE.r12};
% ridx = cellfun(@sum,count) < 200;
% pwcNE(ridx) = [];
% count = {pwcnNE.r12};
% ridx = cellfun(@sum, count) < 200;
% pwcnNE(ridx) = [];
% 
% %get delays
% delay = pwcNE(1).delay;
% ddiff = delay(2)-delay(1);
% dbounds = [delay(1) delay(end)];
% 
% %get number of bins to consider
% maxpd = max([abs([pwcNE.peakdelay]) abs([pwcnNE.peakdelay])]);
% numbins = (dbounds(2) - maxpd) / ddiff;
% 
% %get delay value at pval
% [~, ~, hdNE] = calc_pwc_sharpness_cdf(pwcNE, numbins, pval);
% [~, ~, hdnNE] = calc_pwc_sharpness_cdf(pwcnNE, numbins, pval);
% 

if indivplot == 1 && ~iscell(hdNE)
    error('Indivplot cannot be 1 if you do not have data from more than 1 site')
end

%% Plot individual data
if indivplot == 1
    
    figure;
    c = 1;

    for i = 1:length(hdNE)
        subplot(4,4,c)
        hold on

        boxinput = [hdNE{i};hdnNE{i}];
        group = [zeros(length(hdNE{i}),1);ones(length(hdnNE{i}),1)];

        plotSpread(hdNE{i}, 'distributionColors', [1 .6 .6], 'xValues', 1, 'spreadWidth', 0.5)
        plotSpread(hdnNE{i}, 'distributionColors', [.6 .6 1], 'xValues', 2, 'spreadWidth', 0.5)

        boxplot(boxinput, group, 'notch','on','labels',{sprintf('n=%d',...
            length(hdNE{i})), sprintf('n=%d', length(hdnNE{i}))},...
            'colors','rb', 'symbol', 'b','outliersize',3, 'widths', 0.3);
        
        xlim([0.5 2.5])
        ylim([0 max(hdnNE{i}) + 0.2]);
        
        tickpref;
        
        p(i) = ranksum(hdNE{i}, hdnNE{i});
        
%         if c == 1
%             ylabel(sprintf('Absolute delay at %d^{th} percentile (ms)', 0.5*100))
%         end

%         if c == 16
%             c = 0;
%             figure;
%         end

        c = c + 1;
    end
    print_mfilename(mfilename)

end


% if iscell(hdNE)
%     hdNEmat = cell2mat(hdNE);
% else
%     hdNEmat = hdNE;
% end
% 
% if iscell(hdnNE)
%     hdnNEmat = cell2mat(hdnNE);
% else
%     hdnNEmat = hdnNE;
% end

%statistical test(s)
% p = ranksum(hdNEmat, hdnNEmat);
% fprintf('pval = %.5e\n',p)

hdNEmed = cellfun(@median, hdNE);
hdnNEmed = cellfun(@median, hdnNE);

%plot
figure;
boxinput = [hdNEmed;hdnNEmed];
group = [zeros(length(hdNEmed),1);ones(length(hdnNEmed),1)];

plotSpread(hdNEmed, 'distributionColors', [1 .6 .6], 'xValues', 1, 'spreadWidth', 0.5)
plotSpread(hdnNEmed, 'distributionColors', [.6 .6 1], 'xValues', 2, 'spreadWidth', 0.5)

boxplot(boxinput, group, 'notch','on','labels',{'NE', 'nonNE'},...
    'colors','rb', 'symbol', 'b','outliersize',3, 'widths', 0.2);

ylabel(sprintf('Absolute delay at %d^{th} percentile (ms)', 0.5*100))

print_mfilename(mfilename);

hdNEmat = cell2mat(hdNE);
hdnNEmat = cell2mat(hdnNE);
figure;
boxinput = [hdNEmat; hdnNEmat];
group = [zeros(length(hdNEmat),1); ones(length(hdnNEmat),1)];
plotSpread(hdNEmat, 'distributionColors', [1 .6 .6], 'xValues', 1, 'spreadWidth', 0.5)
plotSpread(hdnNEmat, 'distributionColors', [.6 .6 1], 'xValues', 2, 'spreadWidth', 0.5)

boxplot(boxinput, group, 'notch','on','labels',{'NE', 'nonNE'},...
    'colors','rb', 'symbol', 'b','outliersize',3, 'widths', 0.2);

ylabel(sprintf('Absolute delay at %d^{th} percentile (ms)', 0.5*100))

print_mfilename(mfilename);


end

