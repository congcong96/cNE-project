function Fvals = ne_batch_calc_F1_pvals_spike_subsets(nefiles, subsettype)

Fvals = cell(length(nefiles), 1);

switch subsettype
    case 'i/a'
        
        for i = 1:length(nefiles)

            load(nefiles{i}, 'subset_waveforms')
            if ~exist('subset_waveforms','var')
                continue
            end

            fprintf('\nProcessing F1 values for %s...\n', nefiles{i})
            Fvals{i} = ne_calc_F1_pvals_spike_subsets(subset_waveforms);
            clear('subset_waveforms')

        end 
        
    case 'shared'
        
        for i = 1:length(nefiles)

            load(nefiles{i}, 'shared_waveforms')
            if ~exist('shared_waveforms','var') || isempty(shared_waveforms.projections)
                continue
            end

            fprintf('\nProcessing F1 values for %s...\n', nefiles{i})
            Fvals{i} = ne_calc_F1_pvals_spike_subsets(shared_waveforms);
            clear('shared_waveforms')

        end 
end

    