function pwc = ne_pairwisecorr_function(exp_site_nedata,spktrain,varargin)

if isfield(exp_site_nedata, 'nedata')
    nedata = exp_site_nedata.nedata;
%     df = exp_site_nedata.df;
else
    nedata = exp_site_nedata;
%     df = nedata.df;
end
df = 1;
t = df*0.5;
nbins = size(spktrain,2);
stimdur = t * nbins;

if nargin == 1
    spktrain = nedata.spktrain;
end

try
    mem = nedata.assembly_members;
catch
    mem = nedata.NEmembers;
end


for i = 1:length(mem)
    comb = nchoosek(mem{i},2);
    if isempty(comb)
        pwc{i} = [];
    elseif size(comb,2) == 2
        pwc{i} = pairwisecorr_function(spktrain,stimdur,2000,comb);
    else
        pwc{i} = [];
    end
    
    
end

