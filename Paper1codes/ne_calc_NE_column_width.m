function NEwidth = ne_calc_NE_column_width(exp_site_nedata)

NEmem = exp_site_nedata.nedata.NEmembers;
pos = cell2mat(cellfun(@str2num, exp_site_nedata.nedata.position, 'UniformOutput' ,0));
width = pos(:,1);

NEhorpos = cellfun(@(x) unique(width(x)), NEmem, 'UniformOutput', 0);

NEwidth = cellfun(@(x) max(x) - min(x), NEhorpos);