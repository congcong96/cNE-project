function ne_batch_calc_NE_stamat_with_sta_NEtrain(nefiles)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(nefiles)
    
    fprintf('\nProcessing cNE STAs for %s...\n', nefiles{i})
    load(nefiles{i}, 'exp_site_nedata')
    nedata = exp_site_nedata.nedata;
    NEtrain = nedata.sta_NEtrain;
    
%     dft = 10;
%     dff = 5;
%     rn = exp_site_nedata.stim; % stim as a string
%     rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
%     rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
    nlags = nedata.nlags;
% %     nf = nedata.nf;
%     stimlength = exp_site_nedata.stimlength;
    
    stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);
    exp_site_nedata.nedata.NE_stamat = quick_calc_sta(stimstr.stimulus, NEtrain, nlags);
    
    save(nefiles{i}, 'exp_site_nedata', '-append');
    clear('exp_site_nedata')

end
    


end

