function ne_batch_calc_waveform_silhouette(spkfiles, groupopt)


for i = 1:length(spkfiles)
    
    load(spkfiles{i},'spk','edges')
    basename = regexp(spkfiles{i}, '^\S+(?=(.mat))', 'match', 'once');
    nefile = [basename '-ne-20dft.mat'];
    load(nefile, 'exp_site_nedata');
    
    edges = edges(1:10:end);
    
    
    % Specific to current format of folders and naming conventions. This
    % should be cleaned up at some point.
    try
        if exp_site_nedata.stimlength == 20
            % get rawfolder
            parentrawfolder = 'L:\rawfiles17082018';
            rawfolder = gfn(fullfile(parentrawfolder, sprintf('%s-site%d-*',...
                exp_site_nedata.exp, exp_site_nedata.site)), 1);
            assert(length(rawfolder) == 1);

        elseif exp_site_nedata.stimlength == 30   
            parentrawfolder = 'I:\Data\2017-8-23';
            rawfolder = gfn(fullfile(parentrawfolder, sprintf('site%d_*_frarn1rn16_*_%s',...
                exp_site_nedata.site, exp_site_nedata.exp)), 1);     
            assert(length(rawfolder) == 1);

        else
            try

                parentrawfolder = 'I:\Data\2014-12-15';
                rawfolder = gfn(fullfile(parentrawfolder, sprintf('site%d_*_rn14816rep_%s',...
                    exp_site_nedata.site, exp_site_nedata.exp)), 1); 
                assert(length(rawfolder) == 1);
            catch

                parentrawfolder = 'I:\Data\2015-7-22';
                rawfolder = gfn(fullfile(parentrawfolder, sprintf('site%d_*_frarn14816rep_%s',...
                exp_site_nedata.site, exp_site_nedata.exp)), 1); 
                assert(length(rawfolder) == 1);
            end 

        end
    catch
        warning('Raw folder for %s not found! Skipping...', spkfiles{i});
        continue
        
    end
    
    fprintf('\nProcessing %s spkfile...\n', spkfiles{i});
    
    [waveforms.silhouette, waveforms.projections, waveforms.category] = ne_plot_waveform_pca_clusters(...
        exp_site_nedata, spk, edges, rawfolder{1}, 'groupopt', groupopt, 'plotopt', 0);
    
    switch groupopt
        case 'i/a'
            subset_waveforms_short = waveforms;
            save(nefile,'subset_waveforms_short','-append');
            clear('waveforms', 'subset_waveforms');
        case 'shared'
            shared_waveforms = waveforms;
            save(nefile,'shared_waveforms','-append');
            clear('waveforms','shared_waveforms');
            
    end            
    
end
    