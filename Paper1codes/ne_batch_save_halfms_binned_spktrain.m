function badfiles = ne_batch_save_halfms_binned_spktrain(files, stimfile)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

badfiles = {};

for i = 1:length(files)
    
    fprintf('\nProcessing %s...',files{i})
    
    load(files{i})
    
%     if exist('spktrain','var')
%         fprintf('\nSpiketrain for %s already calculated...\n', files{i})
%         if exist('strf','var')
%             clear('spk','stim_mat','trigger','strf','spktrain','edges','position')
%         else
%             clear('spk','stim_mat','trigger','spktrain','edges','position')
%         end
%     else
    if ~exist('stimfile', 'var')
        
        stim = spk.stim;
%         stimlength = str2double(regexp(spk.stimlength, '^\d{1,3}(?=(min))','match','once'));
        stimlength = 10;


        drive = gcdr;
        subfolder = 'Ripple_Noise\downsampled';
        stimfolder = fullfile(drive,subfolder);

        stimfilename = sprintf('%s-*%dmin_DFt1_DFf5_matrix.mat', stim, stimlength);
        stimfile = fullfile(stimfolder,stimfilename);
        stimfile = gfn(stimfile,1);
        load(stimfile{1}, 'stim_mat');
    else
        load(stimfile{1}, 'stim_mat');
    end
        
%         if isfield(spk,'spk')
%             spk.spk(1).fs = spk.fs;            
        [spktrain, edges] = ne_create_spktrain_from_stim_mat(spk, stim_mat, trigger);
%         else
%             [spktrain, edges, position] = ca_create_spktrain_from_stim_mat(spk, stim_mat, trigger);
%         end
        if isempty(spktrain)
            badfiles = [badfiles files{i}];
            continue
        end

        lastwarn('')

        fprintf('\nSaving %s...\n', files{i})
        if exist('strf','var')
            save(files{i},'spk','trigger','strf','spktrain','edges')
        else
            save(files{i},'spk','trigger','spktrain','edges')
        end
        
        [~, msgid] = lastwarn;
            
        if strcmp(msgid, 'MATLAB:save:sizeTooBigForMATFile')
            if exist('strf','var')
                save(files{i},'spk','trigger','strf','spktrain','edges','-v7.3')
            else
                save(files{i},'spk','trigger','spktrain','edges','-v7.3')
            end
        end
        
%         if exist('strf','var')
%             clear('spk','stim_mat','trigger','strf','spktrain','edges')
%         else
%             clear('spk','stim_mat','trigger','spktrain','edges')
%         end        
        if exist('strf','var')
            clear('spk','stim_mat','trigger','strf','spktrain','edges','stimfile')
        else
            clear('spk','stim_mat','trigger','spktrain','edges','stimfile')
        end
        
%     end
        
end

fprintf('\n')

