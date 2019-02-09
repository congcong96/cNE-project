function [sigpercentage, sigvals] = ne_calculate_ICweight_halves_significance_percentage(corrcell, shuffledcorrcell, prcsig)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Should have used a for loop but now it's quick and unreadable.


sigthresh = cellfun(@(x) prctile(abs(x), prcsig), shuffledcorrcell);
maxcol = cellfun(@(x) max(x,[],1), corrcell, 'UniformOutput', 0);
maxrow = cellfun(@(x) max(x,[],2), corrcell, 'UniformOutput', 0);

collogi = cellfun(@(x,y) x > y, maxcol, num2cell(sigthresh),'UniformOutput',0);
rowlogi = cellfun(@(x,y) x > y, maxrow, num2cell(sigthresh),'UniformOutput', 0);

logivec = [cell2mat(collogi')'; cell2mat(rowlogi)];
sigpercentage = sum(logivec) / length(logivec);

sigcol = cellfun(@(x,y) x(y), maxcol, collogi, 'UniformOutput', 0);
sigrow = cellfun(@(x,y) x(y), maxrow, rowlogi, 'UniformOutput', 0);

sigvals = cell2mat(cellfun(@(x,y) unique([x'; y]), sigcol, sigrow, 'UniformOutput', 0));


end

