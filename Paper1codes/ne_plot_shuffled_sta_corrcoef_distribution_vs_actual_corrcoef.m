function plot_shuffled_sta_corrcoef_distribution_vs_actual_corrcoef...
    (spktrain1, spktrain2, reptrain, inputtype, label, varargin)

narginchk(4, 7)
nlags = 20;
switch inputtype
    case 'spiketrain'
        stim = varargin{1};
        df = varargin{2};


        temp = load(sprintf('%s-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt%d_DFf5_matrix.mat',...
            stim, df));
        stimulus = temp.stim_mat;

        load(sprintf('%s-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt%d_DFf5_param.mat',...
            stim, df), 'faxis');

        acorrval = zeros(size(spktrain1,1),1);
        actsta = zeros(size(spktrain1,1),size(stimulus,1)*nlags);
        modsta = zeros(size(spktrain1,1),size(stimulus,1)*nlags);
        
    case 'sta'
        
        load('rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt1_DFf5_param.mat','faxis');
        
end

ytick = 10:10:length(faxis);
ylab = round(faxis(10:10:length(faxis))/1000); 
figure;
inp = [];

for i = 1:size(spktrain1,1)
    switch inputtype
        case 'spiketrain'
            sta1 = ca_calc_sta_from_stimulus_spktrain(stimulus, spktrain1(i,:), nlags);
            
            sta2 = ca_calc_sta_from_stimulus_spktrain(stimulus, spktrain2(i,:), nlags);
            
%             acorrval(i) = temp(1,2);
        case 'sta'
            sta1 = reshape(spktrain1(i,:), 64, 20);
            sta2 = reshape(spktrain2(i,:), 64, 20);
    
    end
    
    actsta(i,:) = sta1(:);
    modsta(i,:) = sta2(:);
    acorrval(i) = corr(sta1(:), sta2(:));
            
    
    if ~isequal(inp, 'y')
        minmin1 = min(min(sta1));
        maxmax1 = max(max(sta1));
        minmin2 = min(min(sta2));
        maxmax2 = max(max(sta2));
        minmin = min([minmin1 minmin2]);
        maxmax = max([maxmax1 maxmax2]);

        subplot(1,2,1)
        boundary = max([abs(minmin) abs(maxmax)]);
        imagesc(fliplr(sta1));
        xlabel('time before spike (ms)')
        ylabel('frequency (kHz)')
        set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
        set(gca, 'xdir','reverse')
        set(gca,'ydir', 'normal');
        set(gca, 'ytick',ytick,'yticklabel',ylab)
        set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
        set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
%         colormap jet;
        title(sprintf('STA calculated from actual spiketrain of neuron %d', i));

        subplot(1,2,2)
        boundary = max([abs(minmin) abs(maxmax)]);
        imagesc(fliplr(sta2));
        xlabel('time before spike (ms)')
        ylabel('frequency (kHz)')
        set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
        set(gca, 'xdir','reverse')
        set(gca,'ydir', 'normal');
        set(gca, 'ytick',ytick,'yticklabel',ylab)
        set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
        set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
        cmap = cschemes('rdbu', 21);        
        colormap(cmap);
        title(sprintf('STA calculated from modeled spiketrain of neuron %d', i));
%         colorbar;
        
        set(gcf,'Position', [200 200 1500 500])
        
        prompt = 'Done looking at STAs? [y/n]';
        inp = input(prompt,'s');
    end
end

staidx = find(min(abs(acorrval - median(acorrval))) == (abs(acorrval - median(acorrval))), 1);
figure;
sta1 = reshape(actsta(staidx,:), length(faxis), nlags);
sta2 = reshape(modsta(staidx,:), length(faxis), nlags);
% sta3 = reshape(actsta(staidx(2),:), length(faxis), nlags);
% sta4 = reshape(actsta(staidx(2),:), length(faxis), nlags);

minmin1 = min(min(sta1));
maxmax1 = max(max(sta1));
minmin2 = min(min(sta2));
maxmax2 = max(max(sta2));
minmin = min([minmin1 minmin2]);
maxmax = max([maxmax1 maxmax2]);

subplot(1,2,1)
boundary = max([abs(minmin) abs(maxmax)]);
imagesc(fliplr(sta1));
xlabel('time before spike (ms)')
ylabel('frequency (kHz)')
set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
set(gca, 'xdir','reverse')
set(gca,'ydir', 'normal');
set(gca, 'ytick',ytick,'yticklabel',ylab)
set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
% colormap jet;
title(sprintf('STA calculated from actual spiketrain of neuron %d', staidx));

subplot(1,2,2)
boundary = max([abs(minmin) abs(maxmax)]);
imagesc(fliplr(sta2));
xlabel('time before spike (ms)')
ylabel('frequency (kHz)')
set(gca, 'xtick',5:5:20, 'xticklabel',25:25:100)
set(gca, 'xdir','reverse')
set(gca,'ydir', 'normal');
set(gca, 'ytick',ytick,'yticklabel',ylab)
set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
colormap(cmap);
title(sprintf('STA calculated from modeled spiketrain of neuron %d', staidx));

set(gcf,'Position', [200 0 1500 500])

% pause;

% scorrval = zeros(size(spktrain1,1), numiter);

if ~isempty(reptrain)
    
    repcorrval = zeros(size(reptrain,1));
    for j = 1:size(spktrain1,1)
        sta1 = reshape(actsta(j,:), length(faxis), nlags);
        sta2 = ca_calc_sta_from_stimulus_spktrain(stimulus, reptrain(j,:), nlags);
        temp = corrcoef(sta1, sta2);
        repcorrval(j) = temp(1,2);
    end
    
end

% for j = 1:numiter
%     fprintf('\n%d of %d iterations done...', j, numiter)
% 
%     for k = 1:size(spktrain1,1)      
% 
%         surrspkmat = circularly_shuffle_neuronal_spktrain (spktrain1(k,:));
%         surrsta = ca_calc_sta_from_stimulus_spktrain(stimulus, surrspkmat, nlags);
%         temp = corrcoef(actsta(k,:), surrsta);
%         scorrval(k,j) = temp(1,2);
%                
%     end
%     
%     
% end

fprintf('\n')

figure;
hold on
[af, acdf] = ecdf(acorrval(:));
plot(acdf, af, 'r')

set(gcf,'Position', [200 -200 900 700])

xlabel('STA correlation values (x)')
ylabel('F(x)')



if ~isempty(reptrain)
    [sf, scdf] = ecdf(repcorrval(:));
    plot(scdf, sf, 'b')

    legend('Shuffled distribution (control)', 'Model distribution')
    
    [~,p] = kstest2(scdf, acdf);
 
    str = sprintf('P value for KS test:\n%0.4e',p);
    x = xlim;
    x = x(2)-((x(2)-x(1))/5);
    y = ylim;
    y = y(2)-((y(2)-y(1))/4);
    text(x,y, str)

    title(sprintf(['CDF of correlation values between STAs calculated from actual'...
        ' spiketrains and modeled spiketrains (%s)\nand between STAs calculated'...
        ' from actual spiketrains of different presentations of the same stimulus'], label)) 
    
else
    title(sprintf(['CDF of correlation values between STAs calculated from actual'...
        ' spiketrains and modeled spiketrains (%s)'], label)) 
end


hold off

return
