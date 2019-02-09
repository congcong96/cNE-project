function p = ne_compare_spon_evoked_firing_rate(exp_site_nedata)

nedata = exp_site_nedata.nedata;
spktrain = nedata.spktrain;
bidx = nedata.boundaryidx;
NEthresh = nedata.NEthresh;

sspktrain = spktrain(:,1:bidx-1);
espktrain = spktrain(:,bidx:end);

sponrate = sum(sspktrain,2)./(size(sspktrain,2) * 10 / 1000);
evokedrate = sum(espktrain,2)./(size(espktrain,2) * 10 / 1000);

%eliminate outliers temp
idx = sponrate > 10;
sponrate(idx) = [];
evokedrate(idx) = [];

t_act = nedata.total_activities;

NEraster = zeros(size(t_act));

for i = 1:length(NEthresh)
    
    NEraster(i,:) = t_act(i,:) >= NEthresh(i);
    
end

s_act = NEraster(:,1:bidx-1);
e_act = NEraster(:,bidx:end);

NEsponrate = sum(s_act,2)./(size(s_act,2) * 10 / 1000);
NEevokedrate = sum(e_act,2)./(size(e_act,2) * 10 / 1000);

figure; hold on;

for i = 1:size(sponrate,1)
    
    plot(1:2, [sponrate(i) evokedrate(i)],'-o','Color',[0.8,0.8,0.8],'MarkerSize',4)

end

mat = [sponrate evokedrate];
miu = mean(mat);
sem = std(mat)./sqrt(size(mat,1));

errorbar(1:2, miu, sem, 'k', 'LineWidth', 2);
    
xlim([0 3])

tick = 1:2;
set(gca,'xtick', tick, 'xticklabel',{'Spontaneous', 'Evoked'})
ylabel('Neuronal firing rate (Hz)')
tickpref;
print_mfilename(mfilename)

[~, p(1)] = ttest(sponrate, evokedrate);
p(2) = signrank(sponrate, evokedrate);




figure; hold on;

for i = 1:length(NEthresh)
    
    plot(1:2, [NEsponrate(i) NEevokedrate(i)],'-o','Color',[0.8,0.8,0.8],'MarkerSize',4)

end

mat = [NEsponrate NEevokedrate];
miu = mean(mat);
sem = std(mat)./sqrt(size(mat,1));

errorbar(1:2, miu, sem, 'k', 'LineWidth', 2);
    
xlim([0 3])

tick = 1:2;
set(gca,'xtick', tick, 'xticklabel',{'Spontaneous', 'Evoked'})
ylabel('cNE event rate (Hz)')
tickpref;
print_mfilename(mfilename)

[~, p(3)] = ttest(NEsponrate, NEevokedrate);
p(4) = signrank(NEsponrate, NEevokedrate);

