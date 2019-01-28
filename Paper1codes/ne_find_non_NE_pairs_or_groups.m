function [non_NE, thresh, notatall] = ne_find_non_NE_pairs_or_groups(exp_site_nedata, num, numsample, varargin)

% Gets pairs or groups of neurons that are not found in any NE. For groups
% of 3 or larger, as long as 2 of the 3 are found together in a NE, they
% will not be included.

% sets = 4;
% numperset = floor(numsample/4);

p = inputParser;
addParameter(p, 'coincheck', 0, @(x) x == 0 || x == 1)
addParameter(p, 'toinclude', [], @(x) isempty(x) || isscalar(x) || isvector(x))
addParameter(p, 'toincludenum', 1, @isscalar)
addParameter(p, 'flexthresh', 0, @(x) x == 0 || x == 1)
addParameter(p, 'cthresh', 100, @isscalar)
addParameter(p, 'startthresh', 1, @isscalar)
parse(p, varargin{:})
coincheck = p.Results.coincheck;
toinclude = p.Results.toinclude;
toincludenum = p.Results.toincludenum;
flexthresh = p.Results.flexthresh;
cthresh = p.Results.cthresh;
startthresh = p.Results.startthresh;

if ~isempty(toinclude) && toincludenum > length(toinclude)
    error('Number of required neurons included must be smaller or equal to the list of number of neurons to include')
end

if isfield(exp_site_nedata, 'nedata')
    nedata = exp_site_nedata.nedata;
else
    nedata = exp_site_nedata;
end
try
    NEmem = nedata.NEmembers;
catch
    NEmem = nedata.assembly_members;    
end

%get all neurons in NEs
NEmemvec = unique(cell2mat(NEmem));
num_tot = size(nedata.spktrain,1);
tot_neu = 1:num_tot;

%find neurons not in any NE
notatall = find(ismember(tot_neu, NEmemvec) == 0);

unicomb = [];
nsamples = 100000;
thresh = startthresh - 1;
non_NE = [];

c = 1;

% while size(non_NE,1) ~= numsample && thresh < num && (c == 1 || flexthresh == 1)
    
%     c = 1;

    while size(unicomb,1) < numsample && c <= cthresh && (c == 1 || flexthresh == 1)
        
        if flexthresh == 1 || c == 1
            thresh = thresh + 1;
            if flexthresh == 1
                toincludenum = thresh;
            end            
            if thresh > num
                break
            end
            
        end
        

        tot_comb = zeros(nsamples, num);

        % get 100000 random samples of 5 neurons        
        if isempty(toinclude)
            numcomb = nchoosek(num_tot, num);
            
            if numcomb <= 100000
                tot_comb = nchoosek(1:num_tot, num);
            
            else           
            
                for j = 1:nsamples
                    tot_comb(j,:) = sort(randsample(num_tot, num));    
                end
                
            end
            
        else
            try
                numcomb = nchoosek(num_tot - length(toinclude), num - toincludenum) * nchoosek(length(toinclude), toincludenum);
            catch
                if flexthresh == 0
                    non_NE = [];
                    return
                else
                    continue
                end
            end
            
            if numcomb <= 100000
                nontoinclude = setdiff(1:num_tot, toinclude);                             
                sub_tot_comb = nchoosek(nontoinclude, num - toincludenum);
                toincludecomb = nchoosek(toinclude, toincludenum);
                temp = cell(size(toincludecomb,1),1);
                
                for j = 1:size(toincludecomb,1)
                    temp{j} = [repmat(toincludecomb(j,:),size(sub_tot_comb,1),1) sub_tot_comb];
                end
                
                tot_comb = cell2mat(temp);
            
            else            
            
                for j = 1:nsamples
                    if length(toinclude) == toincludenum
                        nottoinclude = setdiff(1:num_tot, toinclude);
                        tot_comb(j,:) = [toinclude sort(randsample(nottoinclude, num - length(toinclude)))];
                    else
                        toincludetemp = randsample(toinclude, toincludenum)';
                        nottoinclude = setdiff(1:num_tot, toincludetemp);
                        tot_comb(j,:) = [toincludetemp sort(randsample(nottoinclude, num - length(toincludetemp)))];
                    end
                end
                
            end
        end
        
        numidx = zeros(size(tot_comb,1), length(NEmem));

        %go through each NE and remove combinations where at least two neurons in
        %each NE is present
        for i = 1:length(NEmem)

            members = NEmem{i};

            numidx(:,i) = sum(ismember(tot_comb, members), 2);

        end
        
        if thresh == 1
            removeidx = sum(numidx >= 2, 2) > 0;
        elseif flexthresh == 0            
            if isempty(toinclude)
                removeidx = sum(numidx > thresh, 2) > 0; %can have more than 1 group of thresh
            else %~isempty(toinclude)
                removeidx = sum(numidx >= 2, 2) > 1 ; %can only have 1 group of thresh (predefined)          
            end      
        else %flexthresh == 1
            removeidx = sum(numidx > thresh, 2) > 0;
        end
        
        tot_comb(removeidx, :) = [];

        %remove groups with 0 coincidence
        if coincheck == 1
            coinratio = ne_calc_coincidence_within_spktrain(nedata.spktrain, tot_comb);
            tot_comb(coinratio == 0,:) = [];
        end

%         unicomb = [unicomb; tot_comb];
        unicomb = unique(tot_comb, 'rows', 'stable');
        c = c+1;
        
        if ~isempty(non_NE)
            [~,commonidx] = intersect(unicomb, non_NE, 'rows');
            unicomb(commonidx,:) = [];      
        end        
        
        if ~isempty(unicomb)            
            
            if size(non_NE,1) + size(unicomb,1) <= numsample
                non_NE = [non_NE; unicomb];
            else
                numreq = numsample - size(non_NE,1);
                non_NE = [non_NE; unicomb(randsample(1:size(unicomb,1),numreq),:)];
            end
        end
        
    end
    
%     if thresh == 1 && flexthresh == 0 && isempty(unicomb)
%         return
%     elseif c <= cthresh && thresh == 1 && size(unicomb,1) > numsample
%         non_NE = unicomb(randsample(1:size(unicomb,1),numsample),:);
%     elseif flexthresh == 1
%         if size(non_NE,1) + size(unicomb,1) < numsample
%             non_NE = [non_NE; unicomb];
%         else
%             non_NE = [non_NE; unicomb(randsample(1:size(unicomb,1),numsample - size(non_NE,1)),:)];
%         end
%     else
%         non_NE = unicomb;
%     end

% end
    
end


