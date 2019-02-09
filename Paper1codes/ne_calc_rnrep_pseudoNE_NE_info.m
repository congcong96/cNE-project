function pseudoNEinfo = ne_calc_rnrep_pseudoNE_NE_info(exp_site_nedata, dtoptim)

% dtoptim in seconds

offset = 200;
fraction = [90 95 97.5 99 100]; % fractions for info calc
nsamples = 10; % number of times to subsample to make num of spks equal in all conditions
niter = 100; % number of times to calc info of frac of sample
spike_count_threshold = 300; % only calc info if each subset has this many spikes
numtrials = 50; %number of trials/repeats
randiter = 50; 

nedata = exp_site_nedata.nedata;
neact = nedata.Activities;
nethresh = nedata.NEthresh;
NEraster = zeros(size(neact));
spktrain = nedata.spktrain;
members = nedata.NEmembers;

for i = 1:length(nethresh)
    NEraster(i,:) = neact(i,:) >= nethresh(i);
end

% initialize subsetinfo struct array
totalNEs = length(members);
pseudoNEinfo(totalNEs).timebin = [];

% stim duration is 4.7s, total duration to consider is 4.7 - offset
% (converted to s)
total_dur = 4.7 - offset/1000;

frac_dur = total_dur * fraction/100;

% calculate pseudo NE rasters
% switch shuffleopt
%     case 'random'

% pNEactthresh = cell(randiter,1);
pNEraster = cell(randiter,1);
pNEact = cell(randiter, 1);

for i = 1:randiter
    pICweights = ne_shuffle_IC_weights(exp_site_nedata, 'random');
    pNEact{i} = assembly_activity(pICweights, spktrain);
%     pNEactthresh{i} = ne_calc_NE_act_thresholds(exp_site_nedata,'circular',20,99,pICweights);
    for j = 1:size(pNEact{i},1)
        pNEraster{i}(j,:) = pNEact{i}(j,:) >= nethresh(j);
%         pNEraster{i}(j,:) = pNEact{i}(j,:) >= pNEactthresh{i}(j);
    end
end

pNEeventcount = cellfun(@(x) sum(x, 2), pNEraster,'UniformOutput',0);
NEeventcount = sum(NEraster,2);
        
c = 1;

for j = 1:length(pseudoNEinfo)

    pseudoNEinfo(c).NE = j;
    minpNEcount = min(cellfun(@(x) x(j), pNEeventcount));
    numcheck = min([NEeventcount(j) minpNEcount]);

    if numcheck > spike_count_threshold
        %% calc NE and neuronal info

        %get minimum number of spikes (95% of the smallest group)
        min_spikes = round(0.95*numcheck);

        %initialize cell array to store all info values                
        irt_NE = cell(nsamples,length(frac_dur));

        for ii = 1:nsamples

            % sub sample spktrain to get the same number of spikes in each group
%             samp_NE = NEraster(j,:);
            samp_NE = sub_sample_spktrain(NEraster(j,:), NEeventcount(j) - min_spikes);

            % group the remaining spiketrains into their respective trials
            trialspk_NE = reshape(samp_NE, [], numtrials);

            %get firing rate of each group
            rt_NE = sum(trialspk_NE,2) ./ numtrials ./ dtoptim;

            %mean firing rate
            rbar_NE = mean(rt_NE);

            for jj = 1:length(frac_dur)

                sampdur = frac_dur(jj); % reduced duration, in seconds
                nvals = floor(sampdur ./ dtoptim); % number of bins in reduced duration

                %calculate information for each fraction
                irt_NE{ii,jj} = calc_rt_subset_info(rt_NE, rbar_NE,...
                    dtoptim, sampdur, nvals, niter);

            end

        end

       %get mean of all info calc for each frac and group
        ifracdurmn_NE = mean(cell2mat(irt_NE));

        %get sd of all info calc for each frac and group
        ifracdursd_NE = std(cell2mat(irt_NE));

        %get extrapolation value for each group
        iextrap_NE = info_extrapolate_from_mn_std(fraction, ifracdurmn_NE);

        %% calculate pseudoNE information

        ifracdurmn_pNE = cell(randiter,1);
        ifracdursd_pNE = cell(randiter,1);
        iextrap_pNE = zeros(randiter,1);

        for kk = 1:randiter

            %initialize cell array to store all info values                
            irt_pNE = cell(nsamples,length(frac_dur));

            for ii = 1:nsamples

                % sub sample spktrain to get the same number of spikes in each group
%                 samp_pNE = pNEraster{kk}(j,:);
                samp_pNE = sub_sample_spktrain(pNEraster{kk}(j,:), pNEeventcount{kk}(j) - min_spikes);

                % group the remaining spiketrains into their respective trials
                trialspk_pNE = reshape(samp_pNE, [], numtrials);

                %get firing rate of each group
                rt_pNE = sum(trialspk_pNE,2) ./ numtrials ./ dtoptim;

                %mean firing rate
                rbar_pNE = mean(rt_pNE);

                for jj = 1:length(frac_dur)

                    sampdur = frac_dur(jj); % reduced duration, in seconds
                    nvals = floor(sampdur ./ dtoptim); % number of bins in reduced duration

                    %calculate information for each fraction
                    irt_pNE{ii,jj} = calc_rt_subset_info(rt_pNE, rbar_pNE,...
                        dtoptim, sampdur, nvals, niter);

                end

            end

            %get mean of all info calc for each frac and group
            ifracdurmn_pNE{kk} = mean(cell2mat(irt_pNE));

            %get sd of all info calc for each frac and group
            ifracdursd_pNE{kk} = std(cell2mat(irt_pNE));

            %get extrapolation value for each group
            iextrap_pNE(kk) = info_extrapolate_from_mn_std(fraction, ifracdurmn_pNE{kk});

        end

        %% save data 


        fprintf('\n');
        fprintf('NE #%d: %.3f bits/spk\n', j, iextrap_NE);
        fprintf('pseudoNE #%d: %.3f bits/spk\n', j, mean(iextrap_pNE));
        fprintf('\n');

        %save data

        pseudoNEinfo(c).timebin = dtoptim;
        pseudoNEinfo(c).fraction = fraction;
        pseudoNEinfo(c).min_spikes = min_spikes;

        pseudoNEinfo(c).NE_info_frac_mn = ifracdurmn_NE;
        pseudoNEinfo(c).NE_info_frac_std = ifracdursd_NE;
        pseudoNEinfo(c).NE_info_extrap = iextrap_NE;

        pseudoNEinfo(c).pseudoNE_info_frac_mn = ifracdurmn_pNE;
        pseudoNEinfo(c).pseudoNE_info_frac_std = ifracdursd_pNE;
        pseudoNEinfo(c).pseudoNE_info_extrap = mean(iextrap_pNE);
        pseudoNEinfo(c).all_pseudoNE_info_extrap = iextrap_pNE;

        c = c+1;

    else
        %if one or more groups did not cross threshold, save empty arrays
        pseudoNEinfo(c).timebin = [];
        pseudoNEinfo(c).fraction = [];
        pseudoNEinfo(c).min_spikes = [];

        pseudoNEinfo(c).NE_info_frac_mn = [];
        pseudoNEinfo(c).NE_info_frac_std = [];
        pseudoNEinfo(c).NE_info_extrap = [];

        pseudoNEinfo(c).pseudoNE_info_frac_mn = [];
        pseudoNEinfo(c).pseudoNE_info_frac_std = [];
        pseudoNEinfo(c).pseudoNE_info_extrap = [];
        pseudoNEinfo(c).all_pseudoNE_info_extrap = [];

        c = c+1;

    end
        
    
end

end

