function ne_plot_ICweights_halves_values(corrmat)

corrvals = cell(length(corrmat), 1);

for i = 1:length(corrmat)
    
    [r, c] = size(corrmat{i});
    
    if r > c
        corrvals{i} = max(corrmat{i}, [], 1);
        corrvals{i} = corrvals{i}(:);
    else %c >= r
        corrvals{i} = max(corrmat{i}, [], 2);
        
    end
end

maxcorrvec = cell2mat(corrvals);

figure;
histogram(maxcorrvec,0:0.05:1,'Normalization','probability');
xlabel('Correlation value')
ylabel('Probability')
x = xlim;
y = ylim;

text(x(2)/4, y(2) - y(2)/4, sprintf('n = %d', length(maxcorrvec)));
tickpref;
print_mfilename(mfilename);


allcorrvec = cell2mat(cellfun(@(x) x(:), corrmat', 'UniformOutput', 0));
minall = round(min(allcorrvec)*10)/10;

figure;
histogram(allcorrvec, minall:0.1:1, 'Normalization', 'probability');
xlabel('Correlation value')
ylabel('Probability')
% 
% yyaxis right
% histogram(maxcorrvec,0:0.1:1,'Normalization','probability');

x = xlim;
y = ylim;

text(x(2)/4, y(2) - y(2)/4, sprintf('n = %d', length(allcorrvec)));
tickpref;
print_mfilename(mfilename);

