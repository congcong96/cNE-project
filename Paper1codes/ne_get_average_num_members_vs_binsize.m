function [nummem, binsize] = ne_get_average_num_members_vs_binsize(fileid, stim, ratioopt)

%Gets average number of members per NE per penetration against binsize
% fileid: Date/time stamp and site number. e.g. 141215_170334-site1

files = gfn([fileid '*' stim '-*dft.mat']);

dft = regexp(files, '(?<=(ne-))\d{1,3}(?=(dft))','match','once');
dft = cellfun(@str2double, dft);

nummem = zeros(length(dft),1);

for i = 1:length(files)
    load(files{i})
    nedata = exp_site_nedata.nedata;
    num = cellfun(@length, nedata.NEmembers);
    nummem(i) = mean(num);
end

[dft, idx] = sort(dft);
binsize = dft * 0.5;
nummem = nummem(idx);

if ratioopt == 1
    nummem = nummem/min(nummem);
end
end

