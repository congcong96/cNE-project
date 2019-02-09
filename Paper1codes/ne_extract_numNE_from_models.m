function ne_extract_numNE_from_models(fileid,stim)

for i = 1:length(fileid)
    
    fprintf('\nProcessing %s...\n', fileid{i})
    
    for j = 1:length(stim)
        
        modelfiles = gfn([fileid{i} '*-' stim{j} '-*NEstats*']);
        base = regexp(modelfiles{1}, '^\S+(?=(-NEstats))','match','once');
        outfile = [base '-numNE'];
        
        try
            expfile = [base '-fs20000-A-spk-strfcmb-ne-20dft.mat'];
            load(expfile)
        catch
            expfile = [base '-fs20000-A-spkcmb-ne-20dft.mat'];
            load(expfile)
        end
        
        numNE.real = size(exp_site_nedata.nedata.Patterns,2);
        
        numNE.preISI = [];
        numNE.preRF = [];
        numNE.prePFRRF = [];
                
        for k = 1:length(modelfiles)
            
            load(modelfiles{k}, 'NEstats')
            preISItemp = [NEstats.preISI_spktrain.numNEs];
            preRFtemp = [NEstats.preRF_spktrain.numNEs];
            prePFRRFtemp = [NEstats.prePFRRF_spktrain.numNEs];
            
%             if k == 1
%                 numNE.real = NEstats.real.numNE;
%             end
            
            numNE.preISI = [numNE.preISI preISItemp];
            numNE.preRF = [numNE.preRF preRFtemp];
            numNE.prePFRRF = [numNE.prePFRRF prePFRRFtemp];
        end
        
        save(outfile,'numNE')
    end
end