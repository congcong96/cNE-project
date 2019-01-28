function stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata, dft)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('dft','var')
    if exp_site_nedata.df > 10
        dft = 10;
    else
        dft = exp_site_nedata.df;
    end
end

curdrive = gcdr;
if dft == 1
    abspath = fullfile(curdrive, 'Ripple_Noise\downsampled');
else
    abspath = fullfile(curdrive, 'Ripple_Noise\downsampled_for_MID');
end

rn = str2double(regexp(exp_site_nedata.stim, '(?<=(rn))\d{1,2}', 'match', 'once'));

if isfield(exp_site_nedata, 'stimlength')
    stimlen = exp_site_nedata.stimlength;
else
    stimlen = 10;
end

dff = 5;

stimstr = ne_get_ripple_noise_stimulus(abspath, rn, dft, dff, stimlen);

end

