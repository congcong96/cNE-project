function [depthmed, depthwidth] = ne_calc_NE_depth_median(exp_site_nedata)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

NEmem = exp_site_nedata.nedata.NEmembers;
pos = cell2mat(exp_site_nedata.nedata.position');
depth = pos(:,2);

depthmed = cellfun(@(x) median(depth(x)), NEmem);

depthwidth = cellfun(@(x) max(depth(x)) - min(depth(x)), NEmem);


