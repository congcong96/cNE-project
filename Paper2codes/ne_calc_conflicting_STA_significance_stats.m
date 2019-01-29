function results = ne_calc_conflicting_STA_significance_stats(NEvet, neuronvet)

sigNE1 = [NEvet.sigNE1];
sigNE2 = [NEvet.sigNE2];
NEinput = [NEvet.input_sig];

results.sigNE1ratio = sum(sigNE1 == NEinput) ./ length(sigNE1);
results.sigNE2ratio = sum(sigNE2 == NEinput) ./ length(sigNE2);
results.ambiguousNEratio = sum(NEinput == 2) ./ length(NEinput);


signeuron1 = [neuronvet.signeuron1];
signeuron2 = [neuronvet.signeuron2];
neuroninput = [neuronvet.input_sig];

results.signeuron1ratio = sum(signeuron1 == neuroninput) ./ length(signeuron1);
results.signeuron2ratio = sum(signeuron2 == neuroninput) ./ length(signeuron2);
results.ambiguousneuronratio = sum(neuroninput == 2) ./ length(neuroninput);