function ne_rn_data_processing(files, DF, roundopt)

% files: cell array of names of files to be processed (typically spkcmb.mat
% or spk-strfcmb.mat)
% DF: vector of DFs (for time) to be processed
% roundopt: if 1, floors downsampled vector to remove the last incomplete
% bin.

% to use this code, binned spiketrains must first be calculated. see 
% ne_batch_save_halfms_binned_spktrain.m

if ~exist('roundopt','var')
    roundopt = 0;
end

drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);

uniquestim = unique(regexp(files, 'rn\d{1,2}-\d{1,2}(?=(min))','match','once'));

if length(uniquestim) == 1
    parts = regexp(uniquestim{1}, '-', 'split');
    
    stimfilename = sprintf('%s-*%smin_DFt1_DFf5_matrix.mat',parts{1}, parts{2});
    stimfile = fullfile(stimfolder,stimfilename);
    stimfile = gfn(stimfile,1);
    load(stimfile{1});
end    


for i = 1:length(files)
    load(files{i});
%     if exist('thresh','var')
%         spk = thresh; %MUA temp
%     end
%     try
%         probetype = spk(1).probetype;
%     catch
%         probetype = 'a1x32-poly3';
%     end

    if ~exist('spktrain','var')
        error('Run ''ne_batch_save_halfms_binned_spktrain.m'' first!')
    end
    
    if isempty(spktrain)
        warning('There are no units in %s. Skipping...', files{i})
        continue
    end
    
    base = regexp(files{i},'^\S+(?=(.mat))','match','once');
    
    if length(uniquestim) ~= 1
        
        stim = regexp(base,'rn\d{1,2}','match','once');
        stimlength = regexp(base, '\d{1,3}(?=(min))', 'match', 'once');

        stimfilename = sprintf('%s-*%smin_DFt1_DFf5_matrix.mat',stim, stimlength);
        stimfile = fullfile(stimfolder,stimfilename);
        stimfile = gfn(stimfile,1);
        load(stimfile{1});
        
    end
    
    position = {spk.spk.position};
    
    for j = 1:length(DF)
        
        DFt = DF(j);
        outfile = sprintf('%s-ne-%ddft.mat',base,DFt);
        outfile = fullfile(pwd,outfile);
        if exist(outfile,'file') == 2
            fprintf('\n%s already exists!\n', outfile)
            continue
        else
            
            nedata = ne_calc_neuronal_ensemble_sta(spktrain, stim_mat, 96000, DFt, edges, position, roundopt);

            save(outfile,'nedata')

            exp_site_nedata = ne_create_exp_site_nedata_file(outfile);

            if isfield(exp_site_nedata.nedata,'CI')
                fprintf('\nCI for %s already calculated\n', outfile);
            else
                CI = ne_calc_ICA_threshold(exp_site_nedata,'circular', 100, 'stdev', 1.5); %threshold currently at 1.5 stdev
                exp_site_nedata.nedata.CI = CI;
            end

            if isfield(exp_site_nedata.nedata,'NEmembers')
                fprintf('\nNEmembers for %s already calculated\n', outfile);
            else
                NEmembers = ne_identify_NEmembers(exp_site_nedata.nedata.Patterns, exp_site_nedata.nedata.CI);
                exp_site_nedata.nedata.NEmembers = NEmembers;

            end

            if isfield(exp_site_nedata.nedata, 'NEthresh')
                fprintf('\ncNE activity thresh for %s already calculated\n', outfile);
            else
                [thresh, alpha] = ne_calc_NE_act_thresholds(exp_site_nedata,'circular', 50, 99:0.1:99.9);
                exp_site_nedata.nedata.NEthresh = thresh;
                exp_site_nedata.nedata.NEthresh_alpha = alpha;
            end

            close all

            save(outfile,'exp_site_nedata');

        
        end
    end
end

end

