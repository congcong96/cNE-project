function ne_plot_constructive_MU_cNE_STA_samples(con_struct, NEfolder)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

assert(length(con_struct) == 1)

if ~exist('NEfolder', 'var')
    NEfolder = 'I:\Cell_Assemblies\MS_NEs_paper2';
end

load(fullfile(NEfolder, con_struct.filename), 'exp_site_nedata')
nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;

% get spike/event trains
MUtrain = logical(sum(nedata.sta_spktrain(con_struct.neurons,:), 1));
NEtrain = logical(nedata.sta_NEtrain(con_struct.NE,:));

% get stimulus
dft = exp_site_nedata.df;
if dft > 10
    dft = 10;
end    
stimlen = exp_site_nedata.stimlength;
dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';

stimstr = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlen);

MUcount = sum(MUtrain);
NEcount = sum(NEtrain);

if NEcount > MUcount
    NEtrain = sub_sample_spktrain(NEtrain, NEcount - MUcount);
else
    MUtrain = sub_sample_spktrain(MUtrain, MUcount - NEcount);
end

sta_NE = calc_single_sta_from_locator_stimulus(NEtrain, stimstr.stimulus, nlags);
sta_MU = calc_single_sta_from_locator_stimulus(MUtrain, stimstr.stimulus, nlags);

boundary = max(abs([sta_NE(:);sta_MU(:)]));

figure('Position', [538 349 860 338]);
cmap = cschemes('rdbu', 100);
colormap(cmap);

subplot(1,2,1)
quick_plot_sta(sta_NE)
title('cNE STRF')

subplot(1,2,2)
quick_plot_sta(sta_MU)
title('Multi-unit STRF')

print_mfilename(mfilename);

end

