function ne_batch_add_trigger_to_spon_evoked_spkfile(spkfiles, trigfolder)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(spkfiles)
    
    load(spkfiles{i})
    
    exp = regexp(spk.exp, '^\d{6}','match','once');
    site = spk.site;
    stim = spk.stim;
    
%     trigfolder = gfn(fullfile(parentfolder, sprintf('site%d_*_%s_*_%s_*', site, stim, exp)), 1);
%     
%     if isempty(trigfolder)
%         exp = num2str(str2double(exp) + 1);
%         % increase date by 1
%         trigfolder = gfn(fullfile(parentfolder, sprintf('site%d_*_%s_*_%s_*', site, stim, exp)), 1);
%         
%         if isempty(trigfolder)
%             warning('\nTrigger file for %s not found!\n', spkfiles{i})
%             continue
%         end
%         
%     end         
    
%     assert(length(trigfolder) == 1)
    triggerfile = gfn(fullfile(trigfolder, sprintf('%s_*-site%d-*-%s-*ADC-00.mat',...
        exp, site, stim)),1);
    load(triggerfile{1}, 'trigger')
    save(spkfiles{i}, 'trigger', '-append')

end

