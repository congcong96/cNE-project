function [tmf, xmf, neuronrtf, NErtf] = ne_calc_ne_neuron_rtf(exp_site_nedata, sigopt)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;

if sigopt
    neuron_sta = nedata.sig_stamat;
else
    neuron_sta = nedata.stamat;
end
nf = nedata.nf;
nlags = nedata.nlags;

dft = 10;
dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';
stimlength = exp_site_nedata.stimlength;

[stimstr] = ne_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlength);
load(stimstr.paramfile, 'taxis','faxis','MaxFM','MaxRD')
taxis = taxis(1:nlags);

modbins = 20;
neuronrtf = zeros(size(neuron_sta,1), (modbins + 1) * (2 * modbins + 1));

for i = 1:size(neuron_sta, 1)
    
    [tmf, xmf, temprtf] = sta2rtf(reshape(neuron_sta(i,:), nf, nlags), taxis, faxis, MaxFM, MaxRD);
    neuronrtf(i,:) = temprtf(:)';
       
end

if sigopt
    NE_sta = nedata.sig_NE_stamat;
else
    NE_sta = nedata.NE_stamat;
end
NErtf = zeros(size(NE_sta,1), (modbins + 1) * (2 * modbins + 1));

for i = 1:size(NE_sta, 1)
    
    [tmf, xmf, temprtf] = sta2rtf(reshape(NE_sta(i,:), nf, nlags), taxis, faxis, MaxFM, MaxRD);
    NErtf(i,:) = temprtf(:)';
    
end

