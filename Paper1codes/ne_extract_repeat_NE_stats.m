function NEstats = ne_extract_repeat_NE_stats(files, repfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

realnumNE = zeros(length(files),1);
realNEsize = zeros(length(files),1);

for i = 1:length(files)
    load(files{i})
    realnumNE(i) = size(exp_site_nedata.nedata.Patterns,2);
    realNEsize(i) = mean(cellfun('length',exp_site_nedata.nedata.NEmembers)); 
end

repfiles = gfn([repfolder '\repNE*']);
repnumNE = zeros(length(repfiles),1);
repNEsize = zeros(length(repfiles),1);
NEidx = cell(length(repfiles),1);

for j = 1:length(repfiles)
    load(repfiles{j})
    repnumNE(j) = repNE.numNEs;
    repNEsize(j) = repNE.NEsize;
    NEidx{j} = repNE.NEidx;
    CI(j,:) = repNE.CI;
end
    

NEstats.realnumNE = realnumNE;
NEstats.realNEsize = realNEsize;
NEstats.repnumNE = repnumNE;
NEstats.repNEsize = repNEsize;
NEstats.repNEidx = NEidx;
NEstats.repCI = CI;