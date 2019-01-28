function ne_rnrep_data_processing(files, totalDF, probetype)

% files: cell array of names of files to be processed (typically spkcmb.mat
% or spk-strfcmb.mat)
% totalDF: vector of totalDF to be processed


for i = 1:length(files)
    load(files{i});
%     probetype = spk(1).probetype;
%     base = regexp(files{i},'^\S+(?=(.mat))','match','once');
%     stim = regexp(base,'rn\d{1,2}rep','match','once');

    if ~exist('spktrain','var')
        error('Run ''ne_batch_save_halfms_binned_spktrain.m'' or other codes that give you binned spike trains first!')
    end
    
    
    for j = 1:length(totalDF)
        basefile = regexp(files{i}, '^\S+(?=(.mat))', 'match', 'once');
        outfile = sprintf('%s-ne-%ddft.mat', basefile, totalDF(j));
        dsspktrain = downsample_spiketrain(spktrain, totalDF(j));
        nedata = ne_data_processing_for_NE_raster(spk, dsspktrain, totalDF(j) * 0.5, 99);
        save(outfile, 'nedata')
        exp_site_nedata = ne_create_exp_site_nedata_file(outfile,probetype);
        save(outfile, 'exp_site_nedata')
        
    end
end
return

        
