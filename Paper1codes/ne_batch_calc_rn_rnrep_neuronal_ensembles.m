function ne_batch_calc_rn_rnrep_neuronal_ensembles(filenames, DFt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);



for i = 1:length(filenames)
    
    longfile = load(filenames{i}, 'spktrain', 'position');
    longspktrain = longfile.spktrain;
    position = longfile.position;
    
    stim = regexp(filenames{i}, '(?<=(db-))\w+(?=(-fs))','match','once');
    repfilename = regexprep(filenames{i}, stim, [stim 'torep']);
    repfilename = regexprep(repfilename, '-strf', '');
    
    repfile = load(fullfile('.\rnrep',repfilename), 'spktrain');
    repspktrain = repfile.spktrain;    
    
    base = regexp(filenames{i},'^\S+(?=(.mat))','match','once');
    stim = regexp(base,'rn\d{1,2}','match','once');
    stimlength = regexp(base, '\d{1,3}(?=(min))', 'match', 'once');
    
    if isempty(stimlength)
        stimlength = num2str(10);
    end
    
    stimfilename = sprintf('%s-*%smin_DFt1_DFf5_matrix.mat',stim, stimlength);
    stimfile = fullfile(stimfolder,stimfilename);
    stimfile = gfn(stimfile,1);
    load(stimfile{1});
       
    base = regexprep(base, stim, [stim stim 'rep']);
    outfile = sprintf('%s-ne-%ddft.mat',base,DFt);
    outfile = fullfile('.\rnrnrep', outfile);
    
    [nedata] = ne_calc_rn_rnrep_neuronal_ensemble_sta(longspktrain, repspktrain, stim_mat, 96000, 20, position);   
    
    save(outfile,'nedata')

    exp_site_nedata = ne_create_exp_site_nedata_file(outfile);
    nedata = exp_site_nedata.nedata;    
    
    CI = ne_calc_ICA_threshold(exp_site_nedata,'circular', 100, 'stdev', 1.5); %threshold currently at 1.5 stdev
    nedata.CI = CI;    
    
    NEmembers = ne_identify_NEmembers(nedata.Patterns, nedata.CI);
    nedata.NEmembers = NEmembers;
    exp_site_nedata.nedata = nedata;
    
    NEthresh = ne_calc_NE_act_thresholds(exp_site_nedata);
    nedata.NEthresh = NEthresh;
    exp_site_nedata.nedata = nedata;
    
    sta_NEtrain = ne_upsample_NEact_using_member_neuron_activity(exp_site_nedata);
    nedata.NEtrain_sta = sta_NEtrain;
    
    sta_stimulus = stim_mat(:,1:10:end);
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(sta_stimulus, sta_NEtrain, nedata.nlags);
    nedata.NEstamat = resp * stim;    
    
    exp_site_nedata.nedata = nedata;
    
    save(outfile, 'exp_site_nedata')
    
    close all
    
end

