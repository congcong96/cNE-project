function NEgroupinfo = ne_calc_ne_randomgroup_sta_info(exp_site_nedata)

pgi = [];
randopt = 'nonNEonly';
memberthresh = 1;
flexthresh = 1;
randiter = 10;
restrictopt = 1;

fraction = [90 95 97.5 99 100];
spike_count_threshold = 300; % only calc info if each subset has this many spikes
nsamples = 10; % number of times to subsample to make num of spks equal in all conditions
niter = 20; % number of times to calc info of frac of sample

nedata = exp_site_nedata.nedata;
nlags = nedata.nlags;
nethresh = nedata.NEthresh;
members = nedata.NEmembers;
NEsize = cellfun('length', members);

if isfield(exp_site_nedata, 'stimlength')
    stimlen = exp_site_nedata.stimlength;
else
    stimlen = 10;
end


if exp_site_nedata.df <= 10
    dft = exp_site_nedata.df;
    spktrain = nedata.spktrain;
else
    dft = 10;
    spktrain = nedata.sta_spktrain;
end

dff = 5;
rn = exp_site_nedata.stim; % stim as a string
rn = str2double(regexp(rn,'(?<=(rn))\d{1,2}','match','once')); % get RN number
rnpath = 'I:\Ripple_Noise\downsampled_for_MID';

stimstr = ca_get_ripple_noise_stimulus(rnpath, rn, dft, dff, stimlen);


% get cNE activity
NEraster = zeros(length(nethresh), size(spktrain,2));
for i = 1:length(nethresh)
    NEraster(i,:) = sum(spktrain(members{i},:), 1);
end

if isempty(pgi)

    thresh = zeros(length(members),1);

    randgrpraster = cell(length(members),1);
    pseudomembers = cell(length(members),1);

    for i = 1:length(NEsize)

        switch randopt
            case 'random'

                for j = 1:randiter

                    pseudomembers{i}(j,:) = randsample(1:size(spktrain,1), NEsize(i));
                    randgrpraster{i}(j,:) = sum(spktrain(pseudomembers{i}(j,:),:));

                end


            case 'nonNEonly'

                if NEsize(i) <= memberthresh
                    continue
                else

                    [pseudomembers{i}, thresh(i)] = ne_find_non_NE_pairs_or_groups(...
                        exp_site_nedata, NEsize(i), randiter, 'toinclude', members{i},...
                        'toincludenum', memberthresh, 'flexthresh', flexthresh, 'cthresh', 1000,...
                        'startthresh', memberthresh);

%                     [pseudomembers{i}, thresh(i)] = ne_find_non_NE_pairs_or_groups(...
%                         exp_site_nedata, NEsize(i), randiter, 'flexthresh', 1, 'cthresh', 1000);

                    for j = 1:size(pseudomembers{i},1)
                        randgrpraster{i}(j,:) = sum(spktrain(pseudomembers{i}(j,:),:));

                    end
                end
        end

    end


    randgrpeventcount = cellfun(@(x) sum(x, 2), randgrpraster,'UniformOutput',0);
    NEeventcount = sum(NEraster,2);


%         remove combinations with less than half or more than double the spike
%         count of the actual NE.
    if restrictopt == 1

        lowerbounds = round(0.5 * NEeventcount);
%         upperbounds = round(2 * NEeventcount);


        for i = 1:length(randgrpeventcount)
            invidx = find(randgrpeventcount{i} < lowerbounds(i)); %| randgrpeventcount{i} > upperbounds(i))

            while ~isempty(invidx)

                switch randopt
                    case 'nonNEonly'

                        temp = [];
                        c = 1;
                        while size(temp,1) < length(invidx)

                            if c == 25
                                break
                            end

                            temp = ne_find_non_NE_pairs_or_groups(exp_site_nedata,...
                                NEsize(i), length(invidx)*2, 'toinclude', members{i},...
                                'toincludenum',memberthresh, 'flexthresh', flexthresh,...
                                'startthresh', memberthresh, 'cthresh', 1000);

                            [~, intidx] = intersect(temp, pseudomembers{i}, 'rows');
                            temp(intidx,:) = [];
                            c = c+1;

                        end

                        if c == 25
                            break
                        end

                        pseudomembers{i}(invidx,:) = temp(randsample(1:size(temp,1), length(invidx)),:);

                        for j = 1:length(invidx) 
                            randgrpraster{i}(invidx(j),:) = sum(spktrain(pseudomembers{i}(invidx(j),:),:));
                        end

                    case 'random'

                        for j = 1:length(invidx)
                            pseudomembers{i}(invidx(j),:) = randsample(1:size(spktrain,1), NEsize(i));
                            randgrpraster{i}(invidx(j),:) = sum(spktrain(pseudomembers{i}(invidx(j),:),:));
                        end
                end
                randgrpeventcount{i} = sum(randgrpraster{i}, 2);
                invidx = find(randgrpeventcount{i} < lowerbounds(i));% | randgrpeventcount{i} > upperbounds(i));
            end

        end

    end

              
else %if ~isempty(pgi)
    
    randgrpraster = cell(length(members),1);
    pseudomembers = cell(length(members),1);
    
    for i = 1:length(pgi)
        
        pseudomembers{i} = pgi(i).pseudomembers;
        
        for j = 1:size(pseudomembers{i},1)
            
            randgrpraster{i}(j,:) = sum(spktrain(pseudomembers{i}(j,:),:));
            
        end
        
    end
    
    randgrpeventcount = cellfun(@(x) sum(x, 2), randgrpraster,'UniformOutput',0);
    NEeventcount = sum(NEraster,2);
    temp = {pgi.member_thresh};
    temp(cellfun('isempty',temp)) = {0};
    thresh = cell2mat(temp);
    
end 

NEgroupinfo(length(NEsize)).NE = [];

for j = 1:length(NEgroupinfo)
    
    NEgroupinfo(j).NE = j;
    minrandgrpcount = min(randgrpeventcount{j});
    numcheck = min([NEeventcount(j) minrandgrpcount]);

    if numcheck > spike_count_threshold && NEsize(j) >= 2 && ~isempty(pseudomembers{j})
        %% calc NE info
        
        min_spikes = round(0.95*numcheck);
        
        %initialize cell array to store all info values                
        NE_ifraction = cell(nsamples,1);

        for i = 1:nsamples

            samp_NE = sub_sample_spktrain(NEraster(j,:), NEeventcount(j) - min_spikes);

            sta_NE = ca_sta_from_locator_stimulus(samp_NE, stimstr.stimulus, nlags);

            [xprior_NE, xposterior_NE] = ne_sta_stimulus_projection(sta_NE, samp_NE, stimstr.stimulus);

            NE_ifraction{i} = ca_subset_info_from_data_fraction2(xprior_NE, xposterior_NE, fraction, niter);

        end

        % Get mean/std of information for the data fractions

        [NE_info_frac_mn, NE_info_frac_std] = ...
            info_from_data_fraction_mean_std(fraction, NE_ifraction);


        % Extrapolate the information values to get the final value for each spike train type
        [NE_info_extrap] = info_extrapolate_from_mn_std(fraction, NE_info_frac_mn);
        
        %% calculate random group information

        randgrp_info_frac_mn = cell(size(pseudomembers{j},1),1);
        randgrp_info_frac_std = cell(size(pseudomembers{j},1),1);
        randgrp_info_extrap = zeros(size(pseudomembers{j},1),1);

        for k = 1:size(pseudomembers{j},1)
            
            %initialize cell array to store all info values                
            randgrp_ifraction = cell(nsamples, 1);
            
            for i = 1:nsamples
                
                samp_randgrp = sub_sample_spktrain(randgrpraster{j}(k,:), randgrpeventcount{j}(k) - min_spikes);
                
                sta_randgrp = ca_sta_from_locator_stimulus(samp_randgrp, stimstr.stimulus, nlags);

                [xprior_randgrp, xposterior_randgrp] = ne_sta_stimulus_projection(sta_randgrp, samp_randgrp, stimstr.stimulus);

                randgrp_ifraction{i} = ca_subset_info_from_data_fraction2(xprior_randgrp, xposterior_randgrp, fraction, niter);                
                
            end

            [randgrp_info_frac_mn{k}, randgrp_info_frac_std{k}] = ...
                info_from_data_fraction_mean_std(fraction, randgrp_ifraction);
            
            [randgrp_info_extrap(k)] = info_extrapolate_from_mn_std(fraction, randgrp_info_frac_mn{k});
            
        end
        
        fprintf('\n');
        fprintf('NE #%d: %.3f bits/spk\n', j, NE_info_extrap);
        fprintf('pseudoNE #%d: %.3f bits/spk\n', j, median(randgrp_info_extrap));
        fprintf('\n');
        
        
        %% save data
        NEgroupinfo(j).fraction = fraction;
        NEgroupinfo(j).min_spikes = min_spikes;
        NEgroupinfo(j).NEmembers = members{j};

        NEgroupinfo(j).NE_info_frac_mn = NE_info_frac_mn;
        NEgroupinfo(j).NE_info_frac_std = NE_info_frac_std;
        NEgroupinfo(j).NE_info_extrap = NE_info_extrap;
        
        NEgroupinfo(j).random_groups = pseudomembers{j};
        NEgroupinfo(j).member_thresh = thresh(j);
        NEgroupinfo(j).randgrp_info_frac_mn = randgrp_info_frac_mn;
        NEgroupinfo(j).randgrp_info_frac_std = randgrp_info_frac_std;
        NEgroupinfo(j).randgrp_info_extrap = median(randgrp_info_extrap);
        NEgroupinfo(j).all_randgrp_info_extrap = randgrp_info_extrap;
        
    else
        
        NEgroupinfo(j).fraction = [];
        NEgroupinfo(j).min_spikes = [];
        NEgroupinfo(j).NEmembers = [];

        NEgroupinfo(j).NE_info_frac_mn = [];
        NEgroupinfo(j).NE_info_frac_std = [];
        NEgroupinfo(j).NE_info_extrap = [];

        NEgroupinfo(j).randgrp_info_frac_mn = [];
        NEgroupinfo(j).randgrp_info_frac_std = [];
        NEgroupinfo(j).randgrp_info_extrap = [];
        
    end
   
end
end