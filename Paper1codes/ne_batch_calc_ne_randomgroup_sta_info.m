function ne_batch_calc_ne_randomgroup_sta_info(files)

for i = 1:length(files)
    
    fprintf('\nProcessing %s...\n', files{i})
    load(files{i})
    
%     if exist('NEgroupinfoshort', 'var')
%         fprintf('\nNEgroupinfo for %s already calculated! Skipping...', files{i})
%         continue
%     end
    
    NEgroupinfoshort = ne_calc_ne_randomgroup_sta_info(exp_site_nedata);
    
%     exp_site_nedata.NEneuroninfo = NEneuroninfo;
    
    save(files{i}, 'exp_site_nedata', 'NEgroupinfoshort')
    clear('exp_site_nedata', 'NEgroupinfoshort')
    
end