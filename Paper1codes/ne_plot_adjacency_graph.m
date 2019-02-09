function [comb, pos] = ne_plot_adjacency_graph(exp_site_nedata, spk, plotopt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

narginchk(2,3)
if nargin == 2
    plotopt = 'all';
end

nedata = exp_site_nedata.nedata;

allpos = [...
    -50 500; 50 500; 
    -50 450; 0 450; 50 450; 
    -50 400; 0 400; 50 400; 
    -50 350; 0 350; 50 350;
    -50 300; 0 300; 50 300;
    -50 250; 0 250; 50 250;
    -50 200; 0 200; 50 200;
    -50 150; 0 150; 50 150;
    -50 100; 0 100; 50 100;
    -50 50; 0 50; 50 50
    -50 0; 0 0; 50 0];

offset = 150;

allposori = [allpos(:,1) exp_site_nedata.depth - offset - allpos(:,2)];

ahoridx = (allposori(:,1) == 0) * 12.5;

atemp = (allposori(:,2) - min(allposori(:,2)))/50;
amultiplier = (((mod(atemp,2)) * 2)- 1)* -1;
ahorshift = round(atemp/2).^2 + round(atemp/2);   
avertidx = 2/3 * amultiplier .* ahorshift;

allpos(:,1) = bsxfun(@plus, allposori(:,1), avertidx);
allpos(:,2) = bsxfun(@plus, allposori(:,2), ahoridx);

pos = reshape([spk.position], length(spk(1).position), length(spk))';
% pos = mat2cell(pos, ones(1,size(pos,1)), 2);
% allposori = mat2cell(allposori, ones(1,size(allposori,1)),2);

[~,shiftidx] = arrayfun(@(x,y) ismember([x y], allposori, 'rows'), pos(:,1),pos(:,2)); 
% horidx = (pos(:,1) == 0) * 12.5;
% 
% temp = (pos(:,2) - min(pos(:,2)))/50;
% multiplier = (((mod(temp,2)) * 2)- 1)* -1;
% horshift = round(temp/2).^2 + round(temp/2);   
% vertidx = 2/3 * multiplier .* horshift;
% 
% pos(:,1) = bsxfun(@plus, pos(:,1), vertidx);
% pos(:,2) = bsxfun(@plus, pos(:,2), horidx);

pos = allpos(shiftidx,:);



NEmem = nedata.NEmembers;

switch plotopt
    case 'eachNE'

        adj = cell(length(NEmem),1);

        for i = 1:length(NEmem)
            comb = nchoosek(NEmem{i}, 2);
            adj{i} = sparse(length(spk),length(spk));
            for j = 1:size(comb,1)
                adj{i}(comb(j,1),comb(j,2)) = 1;
                adj{i}(comb(j,2),comb(j,1)) = 1;
            end

            figure;
            hold on
            gplot(adj{i},pos,'-*')
            scatter(allpos(:,1), allpos(:,2));
            set (gca,'ydir','reverse')
            hold off

        end
        
    case 'all'
        
        adj = sparse(length(spk),length(spk));

        comb = cell2mat(cellfun(@(x) nchoosek(x,2),NEmem, 'UniformOutput', 0));
        comb = unique(comb,'rows');
        for i = 1:size(comb,1)
            adj(comb(i,1),comb(i,2)) = 1;
            adj(comb(i,2),comb(i,1)) = 1;
        end

        figure;
        hold on
        gplot(adj,pos,'-*')
        scatter(allpos(:,1), allpos(:,2));
        xlim([min(allpos(:,1))-10, max(allpos(:,1)) + 10])
        ylim([min(allpos(:,2))-10, max(allpos(:,2)) + 10])
        set (gca,'ydir','reverse')
        title(sprintf('Adjacency graph for %s-site%d-%s',...
            exp_site_nedata.exp, exp_site_nedata.site, exp_site_nedata.stim))
        xlabel('horizontal position w.r.t. center of probe / um')
        ylabel('depth / um')

end

    
end

