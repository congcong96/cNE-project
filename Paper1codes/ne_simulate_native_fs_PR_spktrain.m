function PRspktrain = ne_simulate_native_fs_PR_spktrain(spktrain, stim_mat)

PRspktrain = zeros(size(spktrain));

nf = size(stim_mat,1);
nlags = 100;
fprintf('\nCalculating STAs...\n')
sta = quick_calc_sta(stim_mat, spktrain, nlags, nf);

for i = 1:size(spktrain,1)
    
    fprintf('Processing neuron %d of %d\n', i, size(spktrain,1))
    xprior = ca_calc_projection_from_sta_stimulus_spktrain(reshape(sta(i,:), nf, nlags), stim_mat, spktrain(i,:));
    PRspktrain(i,:) = shufflePVbinsSpikePrediction(xprior, spktrain(i,:), 15, 'random');
    
end