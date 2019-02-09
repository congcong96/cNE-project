function infocell = ne_plot_NEsta_info_vs_binsize(info, plotmode)

% old file format where info struct array has data for all sites

if nargin == 1
    plotmode = 'median';
end

% get sites/stims
fn = fieldnames(info);

% get stims
stim = regexp(fn,'(?<=(site\d{1,2}_))rn\d{1,2}','match','once');

% get dfs
dfs = fieldnames(info.(fn{1}));

% sort fields in increasing dfs
dfsnum = regexp(dfs,'(?<=(^df))\d{1,3}','match','once');
dfsnum = cellfun(@str2double, dfsnum);
[dfsnum,sortidx] = sort(dfsnum);
dfs = dfs(sortidx);
info = structfun(@(x) orderfields(x, dfs), info, 'UniformOutput', 0);

% get binsizes for plotting
binsizes = (dfsnum * 0.5);

boxinput = [];
bsvec = [];
stimcell = [];

for i = 1:length(fn)
    
    data = info.(fn{i});
    
    for j = 1:length(dfs)
        
        temp = [data.(dfs{j}).info_extrap];
        boxinput = [boxinput temp];
        
        % fill stim groups 
        tempstimcell = cell(1, length(temp)); 
        tempstimcell(:) = stim(i);
        stimcell = [stimcell tempstimcell];
        
        % fill binsizes
        bsvec = [bsvec ones(1, length(temp)) * binsizes(j)]; 
        
    end
        
    allinfo(i) = structfun(@(x) [x.info_extrap], data, 'UniformOutput', 0);

end

% plot boxplot without pairing
infocell = squeeze(struct2cell(allinfo))';

figure; hold on

infomean = zeros(size(infocell,2),1);
infosem = zeros(size(infocell,2),1);

for k = 1:size(infocell,2)
    infomat = cell2mat(infocell(:,k)')'; 
    infomed(k) = median(infomat);
    infomean(k) = mean(infomat);
    infosem(k) = std(infomat)/sqrt(length(infomat));
    plotSpread(infomat, 'distributionColors', [.8 .8 .8], 'xValues',...
        binsizes(k), 'spreadWidth', 1.5)
end

switch plotmode
    case 'mean'
        errorbar(binsizes, infomean, infosem, 'k', 'LineWidth', 2);
    case 'median'  
        boxplot(boxinput, bsvec, 'notch','on','labels',cellfun(@num2str,...
            num2cell(binsizes),'UniformOutput', 0), 'colors','b', 'symbol',...
            'b','outliersize',3, 'widths', 1.5, 'position', binsizes)
%         f = fit(binsizes(:),infomed(:),'smoothingspline');       
%         plot(f)
%         plot(binsizes, infomed, 'r')
end

xlabel('Bin size (ms)')
ylabel('Information (bits/event)')
tickpref;

hold off
% [p,~,stats] = kruskalwallis(boxinput, bsvec);
% 
% c = multcompare(stats);
% 

