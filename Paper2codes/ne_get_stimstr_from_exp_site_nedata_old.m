function stimstr = ne_get_stimstr_from_exp_site_nedata_old(exp_site_nedata, dft, rnpath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('dft','var')
    dft = exp_site_nedata.df;
end

if ~exist('rnpath','var')
    rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
end

dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
stimlength = exp_site_nedata.stimlength;

[stimstr] = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlength);

end

