function ne_batch_save_NEgroupsta(files)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(files)
    vars = who('-file',files{i});
    if ismember('NEgroupsta', vars)
        fprintf('\nNEgroupsta for %s already calculated! Skipping...\n', files{i})
        continue
    else
        fprintf('\nProcessing NEgroupsta for %s...\n',files{i})
        load(files{i}, 'exp_site_nedata','NEgroupinfo')
        NEgroupsta = ne_plot_NE_vs_randomgroup_sta(exp_site_nedata, NEgroupinfo, 500, 0);
        save(files{i}, 'NEgroupsta','-append')
        clear('NEgroupsta')
    end
end

