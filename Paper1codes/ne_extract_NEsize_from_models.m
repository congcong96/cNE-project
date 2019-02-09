function ne_extract_NEsize_from_models(fileid,stim)

for i = 1:length(fileid)
    
    for j = 1:length(stim)
        
        modelfiles = gfn([fileid{i} '*-' stim{j} '-*NEstats*']);
        base = regexp(modelfiles{1}, '^\S+(?=(-NEstats))','match','once');
        outfile = [base '-NEsize'];
        try
            expfile = [base '-fs20000-A-spk-strfcmb-ne-20dft.mat'];
            load(expfile)
        catch
            expfile = [base '-fs20000-A-spkcmb-ne-20dft.mat'];
            load(expfile)
        end
        
        
        NEmembers = exp_site_nedata.nedata.NEmembers;
        NEsize.real = mean(cellfun('length', NEmembers));
        
        NEsize.preISI = [];
        NEsize.preRF = [];
        NEsize.prePFRRF = [];
                
        for k = 1:length(modelfiles)
            
            load(modelfiles{k}, 'NEstats')
            preISItemp = [NEstats.preISI_spktrain.NEsize];
            preRFtemp = [NEstats.preRF_spktrain.NEsize];
            prePFRRFtemp = [NEstats.prePFRRF_spktrain.NEsize];
            
%             if k == 1
%                 NEsize.real = NEstats.real.NEsize;
%             end
            
            NEsize.preISI = [NEsize.preISI preISItemp];
            NEsize.preRF = [NEsize.preRF preRFtemp];
            NEsize.prePFRRF = [NEsize.prePFRRF prePFRRFtemp];
        end
        
        save(outfile,'NEsize')
    end
end