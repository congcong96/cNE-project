function ne_calc_surr_spktrain_NEstats(fileid, stim)
clc;

for i = 1:length(fileid)
       
    for j = 1:length(stim)
        
        expfile = gfn([fileid{i} '*' stim{j} '-*20dft.mat']);
        load(expfile{1})
        spktrainlen = size(exp_site_nedata.nedata.spktrain,2);
        files = gfn([fileid{i} '*' stim{j} '-*model*']);
        
        for k = 1:length(files)
            fprintf('Processing %s\n', files{k})
            load(files{k}, 'surrspktrain')
            fn = fieldnames(surrspktrain);
            
            base = regexp(files{k}, '\S+(?=(-model))','match','once');
            outfile = [base '-NEstats-' sprintf('%d.mat', k)];
            
%             outfile = [base '-NEstats.mat'];
            
            if exist(outfile,'file')
                fprintf('%s already exists! Skipping...\n', outfile)
                continue
            end
            
            for ii = 1:length(fn)
                surrcell = surrspktrain.(fn{ii});
                
                for jj = 1:length(surrcell)
                    spktrain = surrcell{jj};
%                     if size(spktrain,2) ~= spktrainlen
%                         diff = size(spktrain,2) - spktrainlen;
%                         spktrain = spktrain(:,1:end-diff);
%                     end
                    
                    ICweights = assembly_patterns(surrcell{jj});
                    
                    if isempty(ICweights)
                        NEstats.(fn{ii})(jj).spktrain = spktrain;
                        NEstats.(fn{ii})(jj).ICwts = [];
                        NEstats.(fn{ii})(jj).numNEs = 0;
                        NEstats.(fn{ii})(jj).CI = [];
                        NEstats.(fn{ii})(jj).NEmembers = [];
                        NEstats.(fn{ii})(jj).NEsize = 0;
                        
                    else
                        
                        CI = ne_calc_ICA_threshold(surrcell{jj}, 'circular', 50, 'stdev', 1.5);

                        NEmembers = cell(size(ICweights,2),1);

                        for kk = 1:size(ICweights,2)

                            idx1 = find(ICweights(:,kk) <= CI(1));
                            idx2 = find(ICweights(:,kk) >= CI(2));
                            NEmembers{kk} = sort([idx1; idx2]);

                        end

                        NEstats.(fn{ii})(jj).spktrain = spktrain;
                        NEstats.(fn{ii})(jj).ICwts = ICweights;
                        NEstats.(fn{ii})(jj).numNEs = size(ICweights,2);
                        NEstats.(fn{ii})(jj).CI = CI;
                        NEstats.(fn{ii})(jj).NEmembers = NEmembers;
                        NEstats.(fn{ii})(jj).NEsize = mean(cellfun('length', NEmembers));
                    end
                   close all 
                end
            end
            save(outfile, 'NEstats', 'surrspktrain')
            clc;
        end
    end
end
