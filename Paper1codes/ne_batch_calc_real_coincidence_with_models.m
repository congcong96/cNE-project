function ne_batch_calc_real_coincidence_with_models(files)
clc;

for i = 1:length(files)
    
    load(files{i})
    
    base = regexp(files{i}, '\S+(?=(-fs))','match','once');

    
    for j = 1:10
        
        outfile = [base sprintf('-model-sameNE-%d.mat',j)];
        
        if exist(outfile, 'file')
            fprintf('%s already exists! Skipping... \n', outfile)
            continue
        end
        
        fprintf('Processing %s...\n', outfile)
        
        if j == 1
            
            [neurons, ratio, NEstats, surrspktrain] = ne_calc_real_coincidence_with_models...
                (exp_site_nedata, 'sameNE', 100, {'preISI','preRF','prePFRRF'},...
                10, 1);
            save(outfile, 'neurons','ratio','NEstats','surrspktrain');

            
        else
            
            if ~exist('neurons', 'var')
                load([base '-model-sameNE-1.mat'], 'neurons')
            end
            
            [neurons, ratio, NEstats, surrspktrain] = ne_calc_real_coincidence_with_models...
                (exp_site_nedata, neurons, 100, {'preISI','preRF','prePFRRF'},...
                10, 1);
            save(outfile, 'neurons','ratio','NEstats','surrspktrain');

            
        end

        clc;
    end
    
    for k = 1:10
        
        outfile = [base sprintf('-model-diffNE-%d.mat',k)];
        
        if exist(outfile, 'file')
            fprintf('%s already exists! Skipping... \n', outfile)
            continue
        end
        
        fprintf('Processing %s...\n', outfile)

        
        if k == 1
            
            [neurons, ratio, NEstats, surrspktrain] = ne_calc_real_coincidence_with_models...
                (exp_site_nedata, 'diffNE', 100, {'preISI','preRF','prePFRRF'},...
                10, 1);
            
            save(outfile, 'neurons','ratio','NEstats','surrspktrain');
            
        else
            
            if ~exist('neurons', 'var')
                load([base '-model-diffNE-1.mat'], 'neurons')
            end
            
             [neurons, ratio, NEstats, surrspktrain] = ne_calc_real_coincidence_with_models...
                (exp_site_nedata, neurons, 100, {'preISI','preRF','prePFRRF'},...
                10, 1);
            
            save(outfile, 'neurons','ratio','NEstats','surrspktrain');

        end
        clc;
    end
    
end
    