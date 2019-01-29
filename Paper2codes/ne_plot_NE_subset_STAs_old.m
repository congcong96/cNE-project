function ne_plot_NE_subset_STAs(exp_site_nedata, NE, sigopt)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;
NEtrain = nedata.sta_NEtrain(NE,:);

if mod(size(NEtrain, 2), 2)
    NEtrain = [NEtrain 0];
    adjustopt = 1;
else
    adjustopt = 0;
end

NEmat = reshape(NEtrain, 2, []);
NEtrainds = sum(NEmat, 1);
members = nedata.NEmembers{NE};

spktrain = nedata.spktrain(members, :);
NEsubset = zeros(size(spktrain,1), size(NEtrain, 2));

for i = 1:size(spktrain, 1)
    
    comptrain = [NEtrainds;spktrain(i,:)];
    templogitrain = comptrain(1,:) == 1 & comptrain(2,:) >= 1;
    
    tempNEtrain = zeros(size(NEmat));
    tempNEtrain(:, templogitrain) = NEmat(:,templogitrain);
    
    NEsubset(i,:) = reshape(tempNEtrain, 1, []);

end

if adjustopt
    NEsubset = NEsubset(:,1:end-1);
end


stimstr = ne_get_stimstr_from_exp_site_nedata(exp_site_nedata);
sta = quick_calc_sta(stimstr.stimulus, NEsubset, nedata.nlags);

if sigopt
    sta_sig = ne_sig_sta_from_stim_obs_resp(sta, NEsubset, stimstr.stimulus, 20, nedata.nlags);
else
    sta_sig = sta;
end

figure;
cmap = flipud(brewermap(1000,'rdbu'));
colormap(cmap);

for i = 1:size(sta_sig, 1)
    subplot(4,4,i)
    quick_plot_sta(sta_sig(i,:)); 
    title(sprintf('cNE subset, neuron #%d\n(number of spikes: %d)',...
        members(i),sum(NEsubset(i,:))));
end