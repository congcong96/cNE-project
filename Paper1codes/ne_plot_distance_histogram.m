function [dist, alldist] = ne_plot_distance_histogram(exp_site_nedata, distopt, plotopt)

if nargin == 1
    plotopt = 1;
    distopt = 'channel';
elseif nargin == 2
    plotopt = 1;
end

NEmem = exp_site_nedata.nedata.NEmembers;

%remove single-member NEs
idx = cellfun(@length, NEmem) == 1;
NEmem(idx) = [];
numneurons = size(exp_site_nedata.nedata.spktrain, 1);


comb = cell2mat(cellfun(@(x) nchoosek(x,2),NEmem, 'UniformOutput', 0));
comb = unique(comb, 'rows');

allcomb = nchoosek(1:numneurons, 2);

pos = cell2mat(exp_site_nedata.nedata.position');

switch distopt
    case 'distance'
        dist = ne_calc_interneuronal_distances(comb, pos, 'distance');
        alldist = ne_calc_interneuronal_distances(allcomb, pos, 'distance');
    case 'channel'
        dist = ne_calc_interneuronal_distances(comb, pos, 'channel');
        alldist = ne_calc_interneuronal_distances(allcomb, pos, 'channel');
end
          
if plotopt == 1
    figure;
    switch distopt
        case 'distance'
            edges = 0:50:hypot(max(pos(:,2)) - min(pos(:,2)), max(pos(:,1)) - min(pos(:,1)))+50;
            histogram(dist,edges)            
        case 'channel'
            histogram(categorical(dist))
    end
    
    title(sprintf('Distance histogram of neurons within NE for %s-site%d-%s',...
    	exp_site_nedata.exp, exp_site_nedata.site, exp_site_nedata.stim))
end

