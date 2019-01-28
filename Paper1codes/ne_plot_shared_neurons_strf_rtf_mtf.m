function res = ne_plot_shared_neurons_strf_rtf_mtf(exp_site_nedata, trigger)

set(0,'DefaultFigureVisible','off')
data = ne_plot_shared_neurons_sta(exp_site_nedata,300,1,1,1);
stamat = data.stamat;
goodcells = data.neurons;
set(0,'DefaultFigureVisible','on')
close all
spkcount = data.spkcount;

for i = 1:length(goodcells)
    
    for j = 1:size(stamat{i},1)
        
        stavec = stamat{i}(j,:);
        n0 = spkcount{i}(j);
        figure;
        res{i}(j) = ca_plot_strf_rtf_mtf(exp_site_cadata, stavec, trigger, n0);
        
      
    end

end

