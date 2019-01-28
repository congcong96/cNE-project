function dg_nedata = ne_ensemble_processing_for_dg_simulated_spktrains(spktrain, niter, CI, dsfactor)

dg_nedata = [];

if ~exist('dsfactor','var')
    dsfactor = 1;
end

if dsfactor > 1
    spktrain = downsample_spiketrain(spktrain, dsfactor);
end

% spktrain = logical(spktrain);

for i = 1:niter
    
    fprintf('\nProcessing iteration %d of %d...\n', i, niter)

    dgspktrain = ne_generate_dg_simulated_spktrains(spktrain);    
    spkmat = downsample_spiketrain(dgspktrain, 20/dsfactor);
    
    ICweights = assembly_patterns(spkmat);
    
    temp.sim_dft = dsfactor;
    temp.spktrain = spkmat;
    temp.patterns = ICweights;
    temp.numNEs = size(ICweights,2);

    
%     if nargin == 4
%         CI = ne_calc_ICA_threshold(temp,'circular', 50, 'stdev', 1.5); %threshold currently at 1.5 stdev
%         temp.CI = CI;
%     else
    temp.CI = CI;
%     end

    NEmembers = ne_identify_NEmembers(temp.patterns, temp.CI);    
    
    temp.NEmembers = NEmembers;
    temp.NEsize = cellfun('length', NEmembers);    
    
    dg_nedata = [dg_nedata temp];
    
    close all

end        
  