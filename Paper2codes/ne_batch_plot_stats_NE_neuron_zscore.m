function ne_batch_plot_stats_NE_neuron_zscore(nefiles, stattype, sig_sta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

NE_zscore = cell(length(nefiles), 1);
neuron_zscore_mu = cell(length(nefiles), 1);
neuron_zscore_sigma = cell(length(nefiles), 1);


for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    
    cs = exp_site_nedata.nedata.sig_classifying_stats;
    members = exp_site_nedata.nedata.NEmembers;
    logicalidx = cellfun('length', members) > 1;

    NE_zscore{i} = cs.(sprintf('NE_%s_zscore',stattype))(logicalidx);
    neuron_zscore_mu{i} = cellfun(@(x) mean(cs.(sprintf('neuron_%s_zscore', stattype))(x)), members(logicalidx));
    neuron_zscore_sigma{i} = cellfun(@(x) std(cs.(sprintf('neuron_%s_zscore', stattype))(x)), members(logicalidx));

end



if ~exist('sig_sta','var')

    figure;
    scatter3(cell2mat(NE_zscore), cell2mat(neuron_zscore_mu), cell2mat(neuron_zscore_sigma), 10, 'filled');
    xlabel(sprintf('NE %s zscore', stattype));
    ylabel(sprintf('Neuron %s zscore mean', stattype));
    zlabel(sprintf('Neuron %s zscore std', stattype));

else
    
    figure; hold on
    NE_zscore = cell2mat(NE_zscore);
    neuron_zscore_mu = cell2mat(neuron_zscore_mu);
    neuron_zscore_sigma = cell2mat(neuron_zscore_sigma);
    
    idx{1} = find([sig_sta.sig_NE] == 0 & [sig_sta.percentage_sig_neurons] > 0.5); %destructive
    idx{2} = find([sig_sta.sig_NE] == 0 & [sig_sta.percentage_sig_neurons] <= 0.5); %stim-independent
    idx{3} = find([sig_sta.sig_NE] == 1 & [sig_sta.percentage_sig_neurons] < 0.5); %constructive
    idx{4} = find([sig_sta.sig_NE] == 1 & [sig_sta.percentage_sig_neurons] >= 0.5); %facilitative
    
    color = eight_color_blind_palette('skyblue','bluegreen','redpurple','vermillion');
    
    for i = 1:length(idx)
        
        scatter3(NE_zscore(idx{i}), neuron_zscore_mu(idx{i}), neuron_zscore_sigma(idx{i}), 10, color(i,:), 'filled');
        
    end
    
    labels = {'Destructive','Stim-independent','Constructive','Facilitative'};
    legend(labels);
    
end

    
    
    
    xlabel(sprintf('NE %s zscore', stattype));
    ylabel(sprintf('Neuron %s zscore mean', stattype));
    zlabel(sprintf('Neuron %s zscore std', stattype));
    
    