function [assemnum, binsize] = ne_get_num_assemblies_vs_binsize(fileid, stim, ratioopt)

%Gets number of cell assemblies per penetration against binsize
% fileid: Date/time stamp and site number. e.g. 141215_170334-site1


files = gfn([fileid '*' stim '-*dft.mat']);

dft = regexp(files, '(?<=(ne-))\d{1,3}(?=(dft))','match','once');
dft = cellfun(@str2double, dft);

assemnum = zeros(length(dft),1);

for i = 1:length(files)
    load(files{i})
    nedata = exp_site_nedata.nedata;
    pat = nedata.Patterns;
    assemnum(i) = size(pat,2);
      
end

[dft, idx] = sort(dft);
binsize = dft * 0.5;
assemnum = assemnum(idx);

if ratioopt == 1
    assemnum = assemnum/min(assemnum);
end