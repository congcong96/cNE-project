function ne_batch_plot_NEsta_info_vs_binsize(files, plotmode)

% new file format: one .mat file for each info site
if nargin == 1
    plotmode = 'median';
end

allinfo = [];
for i = 1:length(files)
    tempinfo = load(files{i});    
    allinfo = [allinfo; structfun(@(x) [x.info_extrap], tempinfo.info, 'UniformOutput', 0)];
end
    
% get dfs
dfs = fieldnames(tempinfo.info);

% sort fields in increasing dfs
dfsnum = regexp(dfs,'(?<=(^df))\d{1,3}','match','once');
dfsnum = cellfun(@str2double, dfsnum);
[dfsnum,sortidx] = sort(dfsnum);
dfs = dfs(sortidx);
allinfo = orderfields(allinfo, dfs);

% get binsizes for plotting
binsizes = (dfsnum * 0.5);

boxinput = [];
bsvec = [];
stimcell = [];

fn = fieldnames(allinfo);

for i = 1:length(fn)
    temp{i} = cell2mat({allinfo.(fn{i})});
    boxinput = [boxinput temp{i}];
    bsvec = [bsvec ones(1, length(temp{i})) * binsizes(i)];     

end

% plot boxplot without pairing
figure; hold on

% infomean = zeros(size(infocell,2),1);
% infosem = zeros(size(infocell,2),1);

for k = 1:length(temp)
%     infomat = cell2mat(infocell(:,k)')'; 
%     infomed(k) = median(infomat);
%     infomean(k) = mean(infomat);
%     infosem(k) = std(infomat)/sqrt(length(infomat));
    plotSpread(temp{k}', 'distributionColors', [.8 .8 .8], 'xValues',...
        binsizes(k)', 'spreadWidth', 1.5)
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
