function nedata = ne_calc_cell_assembly_wo_stim(spk, binsize)

% Get spike train matrix
if isfield(spk, 'position')
    [spktrain, position] = ne_create_spktrain_matrix_wo_stim(spk, binsize);
else
    [spktrain] = ne_create_spktrain_matrix_wo_stim(spk, binsize);
end

% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spktrain);

if exist('position','var')
    ca_plot_cell_assemblies_data(spktrain, Patterns, Activities, position);
    nedata.position = position;
else
    ca_plot_cell_assemblies_data(spktrain, Patterns, Activities);
end


% Assign data to struct for output argument
nedata.spktrain = spktrain;
nedata.Patterns = Patterns;
nedata.Activities = Activities;
nedata.binsize = binsize;
try
    nedata.probe = spk(1).probe;
catch
    warning('No probetype recorded\n')
end