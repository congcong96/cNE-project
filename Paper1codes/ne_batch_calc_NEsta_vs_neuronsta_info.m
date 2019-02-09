function ne_batch_calc_NEsta_vs_neuronsta_info(files)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(files)
    
    filebase = regexp(files{i}, '^\S+(?=(.mat))','match','once');
    outfile = [filebase '-info.mat'];
    load(outfile)
    
%     if ~exist(outfile, 'file')
        clc;
        fprintf('Processing neuron-NE info for %s\n',files{i})
        load(files{i})
        neuronNEinfoproj = ne_calc_NEsta_vs_neuronsta_info(exp_site_nedata, 'project');
%         pseudogroupinfo = ne_calc_NEsta_vs_pseudoNEsta_info(exp_site_nedata);
        save(outfile, 'pseudoNEinfo','neuronNEinfo','pseudogroupinfo', 'neuronNEinfoproj');
        clear('exp_site_nedata','neuronNEinfo', 'pseudoNEinfo','pseudogroupinfo', 'neuronNEinfoproj');
%     else
%         fprintf('%s already processed! Skipping...\n', outfile)
%     end  
    

end

