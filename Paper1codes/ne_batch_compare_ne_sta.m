function [autocorrs, crosscorrs] = ne_batch_compare_ne_sta(nefiles)

for i = 1:length(nefiles)
    fprintf('\nProcessing %s...\n', nefiles{i})
    load(nefiles{i})
    date = regexp(exp_site_nedata.exp, '^\d{6}(?=(_))','match','once');
    site = exp_site_nedata.site;
    stim = exp_site_nedata.stim;
    basefile = regexp(nefiles{i},'\S+(?=(-ne))','match','once');
    spkfile = [basefile '.mat'];
    load(spkfile, 'spktrain');
    spktrain = downsample_spiketrain(spktrain,2);

    [autocorrs.(sprintf('site%d_%s_%s',site,stim,date)), crosscorrs.(sprintf('site%d_%s_%s',site,stim,date))] =...
        ne_compare_ne_sta(exp_site_nedata, spktrain);
    
end