function ne_batch_save_halfms_binned_spktrain_for_spon_evoked_activity(files)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(files)
    
    fprintf('\nProcessing %s...',files{i})
    
    vars = whos('-file',files{i});
    
    if ismember('spktrain', {vars.name})
        fprintf('\nSpiketrain for %s already calculated...', files{i})
    else    
    
        load(files{i}, 'spk', 'trigger')    
        spktimes = {spk.spk.spiketimes};
    
%         stim = spk.stim;
        stim = 'rn1'; %temporarily hardcoded
        stimlength = str2double(regexp(spk.stimlength, '^\d{1,3}(?=(min))','match','once'));

        drive = gcdr;
        subfolder = 'Ripple_Noise\downsampled';
        stimfolder = fullfile(drive,subfolder);

        stimfilename = sprintf('%s-*%dmin_DFt1_DFf5_matrix.mat', stim, stimlength);
        stimfile = fullfile(stimfolder,stimfilename);
        stimfile = gfn(stimfile,1);
        load(stimfile{1});
        
%         if isfield(spk,'spk')
%             spk.spk(1).fs = spk.fs;            
        [evokedspktrain, evokededges, position] = ne_create_spktrain_from_stim_mat(spk, stim_mat, trigger);
%         else
%             [spktrain, edges, position] = ca_create_spktrain_from_stim_mat(spk, stim_mat, trigger);
%         end

        sponend = (trigger(1) - 1)/20000*1000;
        sponedges = 0:0.5:sponend;
        sponspktrain = zeros(length(spk.spk), length(sponedges)-1);
        
        for j = 1:length(spk.spk)
            [sponspktrain(j,:), ~] = histcounts(spktimes{j}, sponedges);
            
        end % (for i)
        
        spktrain = [sponspktrain evokedspktrain];
        edges = [sponedges evokededges];
        
        boundaryidx = length(sponedges);
        

        lastwarn('')

        save(files{i},'spktrain','edges','position','boundaryidx', '-append')
        
        [~, msgid] = lastwarn;
            
        if strcmp(msgid, 'MATLAB:save:sizeTooBigForMATFile')
            
            save(files{i}, 'spk', 'trigger', 'spktrain','edges','position','boundaryidx','-v7.3') %temp without strf
            
        end
        
        clearvars -except files i
        
%         if exist('strf','var')
%             clear('spk','stim_mat','trigger','strf','spktrain','edges','position','boundaryidx')
%         else
%             clear('spk','stim_mat','trigger','spktrain','edges','position','boundaryidx')
%         end        
%         
%         
    end
        
end

fprintf('\n')

