function NEneuronsta = ne_calc_NE_vs_neuron_ptd(exp_site_nedata)

% plotopt: cNE that will be plotted

niter = 10;

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;
members = nedata.NEmembers;
NEraster = nedata.sta_NEtrain;
% sigNE = nedata.sig_NE_sta;
% signeuron = find(nedata.sig_neuron_sta);

% if ~any(sigNE)
%     NEneuronsta = [];
%     return
% end

stimtype = exp_site_nedata.stim;
stimtype = regexp(stimtype, 'rn\d{1,2}','match','once');

spike_count_threshold = 300;
stimlen = exp_site_nedata.stimlength;

if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
    spktrain = nedata.spktrain;
else
    dft = 10;
    spktrain = nedata.sta_spktrain;
end

curdrive = gcdr;

folder = 'Ripple_Noise';
subfolder = 'downsampled_for_MID';

stimmatfile =  gfn(fullfile(curdrive,folder,subfolder,sprintf('%s-*-%dmin_DFt%d_DFf5_matrix.mat',stimtype,stimlen,dft)),1);
load(stimmatfile{1});

% get significant cNEs and their event trains only
% NEraster = NEraster(sigNE, :);
% members = members(sigNE);

% get significant neurons and their spike trains within significant cNEs only
% members = cellfun(@(x) intersect(x, signeuron), members, 'UniformOutput', 0); 

neunespks = cellfun(@(x) logical(spktrain(x,:)), members, 'UniformOutput', 0);


spkcount = cellfun(@(x) sum(x, 2), neunespks, 'UniformOutput', 0);
NEeventcount = sum(NEraster,2);

assert(length(NEeventcount) == length(spkcount))


% initialize info struct array
totalneurons = sum(cellfun(@length, spkcount));

if totalneurons == 0
    NEneuronsta = [];
    return
end

NEneuronsta(totalneurons).neuron = [];

c = 1;

for j = 1:length(spkcount)    
    
    for k = 1:length(spkcount{j})
        
        fprintf('\nCalculating NE #%d vs neuron #%d STA...\n', j, members{j}(k));
        
        NEneuronsta(c).neuron = members{j}(k);
        NEneuronsta(c).NE = j;
        NEneuronsta(c).neuron_count = spkcount{j}(k);
        NEneuronsta(c).NE_count = NEeventcount(j);
        
        min_spikes = min([NEeventcount(j) spkcount{j}(k)]);
        
        if min_spikes > spike_count_threshold
            
            NEneuronsta(c).min_spikes = min_spikes;

            if min_spikes == NEeventcount(j)                
                
                samp_NE = NEraster(j,:);
                samp_neuron = zeros(niter, length(samp_NE));                
                for i = 1:niter
                    samp_neuron(i,:) = sub_sample_spktrain(neunespks{j}(k,:), spkcount{j}(k) - min_spikes);
                end
                
                NEneuronsta(c).sta_NE = calc_single_sta_from_locator_stimulus(samp_NE, stim_mat, nlags, 1);
                NEneuronsta(c).sta_neuron = quick_calc_sta(stim_mat, samp_neuron, nlags);
                NEneuronsta(c).ptd_NE = max(NEneuronsta(c).sta_NE) - min(NEneuronsta(c).sta_NE);
                NEneuronsta(c).ptd_neuron = max(NEneuronsta(c).sta_neuron, [], 2) - min(NEneuronsta(c).sta_neuron, [], 2);
                
            else
                samp_neuron = neunespks{j}(k,:);
                samp_NE = zeros(niter, length(samp_neuron)); 
                for i = 1:niter
                    samp_NE(i,:) = sub_sample_spktrain(NEraster(j,:), NEeventcount(j) - min_spikes);
                end
                
                NEneuronsta(c).sta_NE = quick_calc_sta(stim_mat, samp_NE, nlags);
                NEneuronsta(c).sta_neuron = calc_single_sta_from_locator_stimulus(samp_neuron, stim_mat, nlags, 1);
                NEneuronsta(c).ptd_NE = max(NEneuronsta(c).sta_NE, [], 2) - min(NEneuronsta(c).sta_NE, [], 2);
                NEneuronsta(c).ptd_neuron = max(NEneuronsta(c).sta_neuron) - min(NEneuronsta(c).sta_neuron);
            
            end            

            c = c+1;
        else
            NEneuronsta(c).min_spikes = min_spikes;
            c = c+1;
            continue
        end
    end
end

end

