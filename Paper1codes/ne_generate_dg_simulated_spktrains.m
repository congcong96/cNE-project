function dgspktrain = ne_generate_dg_simulated_spktrains(spktrain)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% spktrain = logical(spktrain);
sigma = cov(spktrain');
mu = mean(spktrain,2);
nsamples = size(spktrain,2);

dgspktrain = (sampleDichGauss01(mu,sigma,nsamples)');

end

