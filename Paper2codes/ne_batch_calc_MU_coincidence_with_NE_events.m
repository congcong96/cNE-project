function [NEcoincidence, PSTH] = ne_batch_calc_MU_coincidence_with_NE_events(nefiles)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

MUfolder = 'I:\Multi-unit';
prevstim = '';
prevstimlen = NaN;

NEcoincidence = cell(length(nefiles), 1);
PSTH = cell(length(nefiles), 1);

for i = 1:length(nefiles)
    
    fprintf('\nProcessing %s (file %d of %d)...\n', nefiles{i}, i, length(nefiles))
    
    load(nefiles{i}, 'exp_site_nedata')
    
    basefile = regexp(nefiles{i}, '\S+(?=-spk-filt)','match','once');
    threshfile = fullfile(MUfolder, [basefile '-thresh.mat']);
    try
        load(threshfile, 'thresh', 'trigger')
    catch
        warning('\nMU file for %s not found! Skipping...\n', nefiles{i})
        continue
    end
    
    if ~(strcmp(prevstim, exp_site_nedata.stim) && prevstimlen == exp_site_nedata.stimlength)   
        stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata, 1);
    end    
    
    [NEcoincidence{i}, PSTH{i}] = ne_calc_MU_coincidence_with_NE_events(exp_site_nedata, thresh, trigger, stimstr.stimulus);
    
    prevstim = exp_site_nedata.stim;
    prevstimlen = exp_site_nedata.stimlength;
    
end
    