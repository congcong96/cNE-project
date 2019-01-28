function ne_batch_calc_NE_for_spon_evoked_activity(files, DF)

for i = 1:length(files)    
    
    basename = regexp(files{i}, '^\S+(?=(.mat))','match','once');    
    outfile = [basename sprintf('-ne-%ddft.mat', DF)];
    
%     if ~exist(outfile,'file')
%         
%         vars = who('-file', outfile);
% 
% 
%         if ismember('nedata', vars)
% 
%             fprintf('\nnedata for %s already calculated. Skipping...', outfile)
%             load(outfile, 'nedata')
% 
%         else


            load(files{i}, 'spktrain','edges','boundaryidx','position')

            fprintf('\nProcessing cNEs for %s...', outfile)

            nedata = ne_calc_NE_for_spon_evoked_activity(spktrain, DF, edges, boundaryidx, position);

            save(outfile, 'nedata');

            close all

%         end
%     end
%     
    exp_site_nedata = ne_create_exp_site_nedata_file(outfile);
    nedata = exp_site_nedata.nedata;

    CI_all = ne_calc_ICA_threshold(nedata.spktrain,'circular', 100, 'stdev', 1.5); %threshold currently at 1.5 stdev
    nedata.CI_all = CI_all;
    
    CI_spon = ne_calc_ICA_threshold(nedata.spktrain(:,1:nedata.boundaryidx-1), 'circular', 100, 'stdev', 1.5);
    nedata.CI_spon = CI_spon;
    
    CI_evoked = ne_calc_ICA_threshold(nedata.spktrain(:,nedata.boundaryidx:end), 'circular', 100, 'stdev', 1.5);
    nedata.CI_evoked = CI_evoked;
    
    close all
    
    NEmembers_all = ne_identify_NEmembers(nedata.total_patterns, CI_all);
    nedata.NEmembers_all = NEmembers_all;
    
    NEmembers_spon = ne_identify_NEmembers(nedata.spon_patterns, CI_spon);
    nedata.NEmembers_spon = NEmembers_spon;
    
    NEmembers_evoked = ne_identify_NEmembers(nedata.evoked_patterns, CI_evoked);
    nedata.NEmembers_evoked = NEmembers_evoked;
  
    exp_site_nedata.nedata = nedata;
    
    save(outfile, 'exp_site_nedata');
    
    fprintf('\n')

    
end