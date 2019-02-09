function [NEmemdepth, neudepth] = ne_calc_NEmembers_depth(exp_site_nedata)

NEmem = unique(cell2mat(exp_site_nedata.nedata.NEmembers));
pos = cell2mat(exp_site_nedata.nedata.position');
neudepth = pos(:,2);

NEmemdepth = neudepth(NEmem);
