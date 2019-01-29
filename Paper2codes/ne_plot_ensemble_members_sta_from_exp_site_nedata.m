function ne_plot_ensemble_members_sta_from_exp_site_nedata(exp_site_nedata, sigopt)

% 
% if ~exist('ylabelopt', 'var')
%     ylabelopt = 'freq';
% end
% 
% if ~exist('ylimopt', 'var')
%     ylimopt = [];
% end

nedata = exp_site_nedata.nedata;
NEmembers = nedata.NEmembers;
NEact = nedata.Activities;
spktrain = nedata.sta_spktrain;
NEthresh = nedata.NEthresh;
ICwt = nedata.Patterns;
CI = nedata.CI;

filename = sprintf('%s-site%d-%s', exp_site_nedata.exp, exp_site_nedata.site, exp_site_nedata.stim);


if ~isfield(exp_site_nedata, 'stimlength')
    stimlength = 10;
else
    stimlength = exp_site_nedata.stimlength;
end

if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
else
    dft = 10;
end

dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
nlags = nedata.nlags;

[stimstr] = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlength);

if ~isfield(nedata, 'NE_stamat')
    NEtrain = ne_upsample_NEact_using_member_neuron_activity(NEact, NEmembers, spktrain, NEthresh);
    NEsta = quick_calc_sta(stimstr.stimulus, NEtrain, nlags);
else
    if sigopt
        NEsta = nedata.sig_NE_stamat;
    else
        NEsta = nedata.NE_stamat;
    end
end

load(stimstr.paramfile, 'faxis')
ytick = 20:20:length(faxis);
ylab = round(faxis(20:20:length(faxis))/1000);


if ~isfield(nedata,'NEmembers')
    error('Please calculate threshold and identify members of cell assembly first');
end

if sigopt
    stamat = nedata.sig_stamat;     
else
    stamat = nedata.stamat;
end


for i = 1:size(NEsta,1)
    
    ne_plot_ensemble_members_sta(ICwt(:,i), CI, NEmembers{i}, stamat, ...
        NEsta(i,:), i, 'separate', 'ytick', ytick, 'ylab', ylab, ...
        'filename', filename);
    
end