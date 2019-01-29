function ne_batch_calc_ne_neuron_sta_info(files)

for i = 1:length(files)
    
    fprintf('\nProcessing %s...\n', files{i})
    load(files{i})
    
%     if exist('NEneuroninfo', 'var')
%         
%         fprintf('\nNEneuroninfo for %s already calculated! Skipping...', files{i})
%         clear('NEneuroninfo')
%         continue
% %         
% %     elseif isfield(exp_site_nedata, 'NEneuroninfo')
% %         
% %         NEneuroninfo = exp_site_nedata.NEneuroninfo;
% %         exp_site_nedata = rmfield(exp_site_nedata, 'NEneuroninfo');
% %         save(files{i}, 'exp_site_nedata','NEneuroninfo','-append')
% %         
%     else
       
    NEneuroninfo = ne_calc_ne_neuron_sta_info(exp_site_nedata, 1, NEneuroninfo);
    save(files{i}, 'NEneuroninfo', '-append')
        
%     end
    
    clear('exp_site_nedata', 'NEneuroninfo')
    
end