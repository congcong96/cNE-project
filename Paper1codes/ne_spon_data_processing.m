function ne_spon_data_processing (files, binsize)

if ischar(files)
    files = {files};
end


for i = 1:length(files)
    
    load(files{i})
%     [spktrain, position, edges] = ca_create_spktrain_matrix_wo_stim(spk, 0.5);
%     stimdur = size(spktrain,2) * 0.5;
%     allpwc = pairwisecorr_function(spktrain,stimdur);
%     
%     save(files{i}, 'spktrain', 'spk', 'position', 'edges','allpwc')

    if isempty(spk.spk)
        continue
    end
    
    base = regexp(files{i},'^\S+(?=(.mat))','match','once');
    
    oribinsize = 0.5;
    dft = binsize/oribinsize;

    outfile = sprintf('%s-ne-%ddft.mat',base,dft);
    
    if exist(outfile, 'file')
       continue
    end
    
    fprintf('\nProcessing %s...\n', outfile)
    
    [nedata] = ne_calc_cell_assembly_wo_stim(spk, binsize);    

    save(outfile,'nedata')

    exp_site_nedata = ne_create_exp_site_nedata_file(outfile);
    nedata = exp_site_nedata.nedata;

        if isfield(nedata,'CI')
            fprintf('\nCI for %s already calculated\n', outfile);
        else
            CI = ne_calc_ICA_threshold(exp_site_nedata,'circular');
            nedata.CI = CI;
            exp_site_nedata.nedata = nedata;
        end

        if isfield(nedata,'NEmembers')
            fprintf('\naNEmembers for %s already calculated\n', outfile);
        else
            NEmembers = ne_identify_NEmembers(exp_site_nedata.nedata.Patterns, exp_site_nedata.nedata.CI);
            nedata.NEmembers = NEmembers;
            exp_site_nedata.nedata = nedata;

        end

%         if isfield(nedata,'pwc')
%             fprintf('\nPairwise correlations for %s already calculated\n', outfile);
%         else
% 
% %             [spktrain] = ca_create_spktrain_matrix_wo_stim(spk, 0.5);
%             pwc = ne_pairwisecorr_function(exp_site_nedata,spktrain, allpwc);
%             nedata.pwc = pwc;
%             exp_site_nedata.nedata = nedata;
% 
% 
%         end
    close all

    save(outfile,'exp_site_nedata');

end



