function ne_plot_coinratio_explained_vs_num_neurons(ratio, mtype)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('mtype','var')
    mtype = 'median';
end

figure;
hold on
fn = fieldnames(ratio);

% find real ratios
rridx = ~cellfun('isempty', strfind(fn, 'real'));
realratio = ratio.(fn{rridx});

otherratios = rmfield(ratio, 'real');
labels = fieldnames(otherratios);
colors = distinguishable_colors(length(labels), {'w','k'});
% colors([4 5],:) = colors([5 4],:);

for i = 1:length(labels)
    ratiocell = ratio.(labels{i});
    normcell = cellfun(@(x,y) bsxfun(@rdivide, x, y), ratiocell, realratio, 'UniformOutput', 0);
    
    switch mtype
        case 'median'
            yvals = cellfun(@(x) median(x, 'omitnan'), normcell);
        case 'mean'
            yvals = cellfun(@(x) mean(x, 'omitnan'), normcell);
    end
    
    xvals = 2:length(yvals)+1;
    
%     upper = cellfun(@(x) prctile(x, 97.5), normcell);
%     lower = cellfun(@(x) prctile(x, 2.5), normcell);
    
    plot(xvals,yvals,'Color',colors(i,:))
%     plot(xvals,upper,'Color',colors(i,:),'LineStyle','--')
%     plot(xvals,lower,'Color',colors(i,:),'LineStyle','--')

    

    xlabel('Set size')
    ylabel('Median normalized coincidence ratio')


end

legend(labels)

line([min(xvals) max(xvals)],[1 1],'Color','k','LineStyle','--')
ylim([0 1.1])
tickpref;
print_mfilename(mfilename);