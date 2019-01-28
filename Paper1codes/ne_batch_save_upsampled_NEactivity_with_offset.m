function ne_batch_save_upsampled_NEactivity_with_offset(NEfiles)

rn16stimmat = load('rn16-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt1_DFf5_matrix.mat');
rn1stimmat = load('rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt1_DFf5_matrix.mat');


for i = 1:length(NEfiles)
    clc;
    fprintf('Processing %s...\n', NEfiles{i})
    load(NEfiles{i})
    df = exp_site_nedata.df;
    stim = exp_site_nedata.stim;

    basefilename = regexp(NEfiles{i}, '^\S+(?=(-spk))','match','once');
    spkfile = [basefilename '-spk-strfcmb.mat'];
    
    switch stim
        case 'rn1'
            [NEraster, usNE, NEact, NEthresh] = ne_upsample_NEactivity_with_offset(spkfile, df, 10, rn1stimmat.stim_mat);
        case 'rn16'
            [NEraster, usNE, NEact, NEthresh] = ne_upsample_NEactivity_with_offset(spkfile, df, 10, rn16stimmat.stim_mat);
    end
    
    
     % make NE raster the same length as neuronal spiketrain
    nedata = exp_site_nedata.nedata;
    spktrain = nedata.sta_spktrain;
    sizediff = size(NEraster,2) - size(spktrain,2);
    NEraster = NEraster(:,1:end-sizediff);
    NEact = NEact(:,1:end-sizediff);
    
    % match NEs
    corrmat = abs(corr(usNE(1).Patterns, nedata.Patterns));
    [row, col] = find(corrmat > 0.99);
    
%     temppattern = usNE(1).Patterns(:, row);
    
    MID_NEtrain = cell(size(nedata.Patterns,2), 1);
    sta_NEtrain = cell(size(nedata.Patterns,2), 1);
    sta_NEthresh = cell(size(nedata.Patterns,2), 1); 
    temppat = cell(1, size(nedata.Patterns,2));

    for j = 1:length(col)
        temppat{col(j)} = usNE(1).Patterns(:,row(j));
        MID_NEtrain{col(j)} = NEraster(row(j),:);
        sta_NEtrain{col(j)} = NEact(row(j),:);
        sta_NEthresh{col(j)} = NEthresh(row(j));
    end
    
    exp_site_nedata.nedata.pattern_corrmat = corr(cell2mat(temppat), nedata.Patterns);
    exp_site_nedata.nedata.MID_NEtrain = MID_NEtrain;
    exp_site_nedata.nedata.sta_NEtrain = sta_NEtrain;
    exp_site_nedata.nedata.sta_NEthresh = sta_NEthresh;   
    
    save(NEfiles{i}, 'exp_site_nedata')
end
