function varargout = ne_plot_NE_subset_STAs(exp_site_nedata, NE, sigopt, type, plotopt)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('plotopt','var')
    plotopt = 0;
end

nedata = exp_site_nedata.nedata;
NEtrain = nedata.sta_NEtrain(NE,:);

members = nedata.NEmembers{NE};
spktrain = nedata.spktrain(members, :);

NEsubset = ne_get_NEsubsets_from_member_spiketrain(NEtrain, spktrain, type);

stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);
sta = quick_calc_sta(stimstr.stimulus, NEsubset, nedata.nlags);

if sigopt
    sta_sig = ne_sig_sta_from_stim_obs_resp(sta, NEsubset, stimstr.stimulus, 20, nedata.nlags);
else
    sta_sig = sta;
end

cmap = flipud(brewermap(1000,'rdbu'));

switch type
    case 'individual_members'
        
        if plotopt
            figure;
            colormap(cmap);

            for i = 1:size(sta_sig, 1)
                subplot(4,4,i)
                quick_plot_sta(sta_sig(i,:)); 
                title(sprintf('cNE subset, neuron #%d\n(number of spikes: %d)',...
                    members(i),sum(NEsubset(i,:))));
            end
        end
        
%         print_mfilename(mfilename);
        
    case 'with/without'
        
        if plotopt
        
            figure('Position', [435 561 1201 315]);
            colormap(cmap);

            subplot(131)
            if sigopt
                quick_plot_sta(exp_site_nedata.nedata.sig_NE_stamat(NE,:));
            else
                quick_plot_sta(exp_site_nedata.nedata.NE_stamat(NE,:));
            end
            title(sprintf('cNE all events\n(number of spikes: %d)', sum(NEsubset(:))))


            subplot(132)
            quick_plot_sta(sta_sig(1,:));
            title(sprintf('cNE subset, with any\n(number of spikes: %d)',...
                sum(NEsubset(1,:))))

            subplot(133)
            quick_plot_sta(sta_sig(2,:));
            title(sprintf('cNE subset, without\n(number of spikes: %d)',...
                sum(NEsubset(2,:))))
        end
        
        if sigopt
            corrval(1) = corr(exp_site_nedata.nedata.sig_NE_stamat(NE,:)', sta_sig(1,:)');
            corrval(2) = corr(exp_site_nedata.nedata.sig_NE_stamat(NE,:)', sta_sig(2,:)');
        else
            corrval(1) = corr(exp_site_nedata.nedata.NE_stamat(NE,:)', sta_sig(1,:)');
            corrval(2) = corr(exp_site_nedata.nedata.NE_stamat(NE,:)', sta_sig(2,:)');
        end
        
        varargout{1} = corrval;

%         print_mfilename(mfilename);
end