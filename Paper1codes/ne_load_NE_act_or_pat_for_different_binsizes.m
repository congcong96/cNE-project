function NEdata = ne_load_NE_act_or_pat_for_different_binsizes(id,stim,df,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

files = gfn(sprintf('%s-*-%s-*dft.mat',id, stim));

if ~isempty(df)
    idx = [];
    for i = 1:length(df)
        temp = regexp(files, sprintf('(?<=(ne-))%ddft.mat',df(i)), 'match', 'once');
        tempidx = find(~cellfun('isempty', temp));
        idx = [idx tempidx];
    end
    
    files = files(idx);
end
   
    
for j = 1:length(files)

    load(files{j})
    df = exp_site_nedata.df;
    switch opt
        case 'act'
            NEdata.(sprintf('df%d', df)) = exp_site_nedata.nedata.Activities;
        case 'pat'
            NEdata.(sprintf('df%d', df)) = exp_site_nedata.nedata.Patterns;
    end
end


