function ne_batch_NE_processing_for_dg_simulated_spktrains(nefiles, dft)

numbatch = 20;
numiter = 10;

for i = 1:length(nefiles)
    
    clc;
    fprintf('Processing %s...\n', nefiles{i})
    
    for j = 1:numbatch
        
        fprintf('\nBatch %d of %d...\n', j, numbatch);

        load(nefiles{i}, 'exp_site_nedata');

        base = regexp(nefiles{i},'^\S+(?=(-ne))','match','once');
        outfile = sprintf('%s-dsdgsim-fixedCI-%d.mat', base, j);
        
        if exist(outfile,'file')
            fprintf('\n%s already exists! Skipping...', outfile)
            continue
        end
        
        spkfile = [base '.mat'];
        load(spkfile, 'spktrain');

        outfile = sprintf('%s-dsdgsim-fixedCI-%d.mat', base, j);
        outfile = fullfile(pwd,outfile);        
        
        if exist(outfile, 'file')
            continue
        end
        
%         spktrain = exp_site_nedata.nedata.sta_spktrain;
        CI = exp_site_nedata.nedata.CI;        

        dg_nedata = ne_ensemble_processing_for_dg_simulated_spktrains(spktrain, numiter, CI, dft);


        save(outfile, 'dg_nedata');
    end
end

