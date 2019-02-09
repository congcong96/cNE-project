function [ratio, varargout] = ne_calc_coincidence_from_spktrains(exp_site_nedata, comb, othertrains, numfiles, mtype)

realspktrain = exp_site_nedata.nedata.spktrain;

for i = 1:length(comb)
    fprintf('\nProcessing sets of %d...\n', i+1)
    realratio{i} = ne_calc_coincidence_within_spktrain(realspktrain, comb{i});
    [repratio{i}, allrepratio{i}] = ne_calc_coincidence_across_spktrains(realspktrain, othertrains, comb{i}, 'allcomb', 1000, mtype);
    
end

cd('I:\Cell_Assemblies\PaperNEsnewthresh\stimreps\modelled_spiketrains')

preISIfiles = gfn('*shuffle*');
preISIfiles = preISIfiles(1:numfiles);
[preISIratio, allISIratio] = ne_calc_coinratio_from_simulated_spktrains(preISIfiles, 'preISI', comb, mtype);

preRFfiles = gfn('*-PR*');
preRFfiles = preRFfiles(1:numfiles);
[preRFratio, allRFratio] = ne_calc_coinratio_from_simulated_spktrains(preRFfiles, 'preRF', comb, mtype);

% prePFRRFfiles = gfn('*-PPR*');
% prePFRRFfiles = prePFRRFfiles(1:numfiles);
% [prePFRRFratio, allPFRRFratio] = ne_calc_coinratio_from_simulated_spktrains(prePFRRFfiles, 'prePFRRF', comb, mtype);

DGfiles = gfn('*-ds4DG*');
DGfiles = DGfiles(1:numfiles);
[DGratio, allDGratio] = ne_calc_coinratio_from_simulated_spktrains(DGfiles, 'DG', comb, mtype);

ratio.real = realratio;
ratio.repeat = repratio;
ratio.shuffle = preISIratio';
ratio.PR = preRFratio';
% ratio.PPR = prePFRRFratio';
ratio.DG = DGratio';

allratio.repeat = allrepratio;
allratio.shuffle = allISIratio;
allratio.PR = allRFratio;
% allratio.PPR = allPFRRFratio;
allratio.DG = allDGratio;

varargout{1} = allratio;