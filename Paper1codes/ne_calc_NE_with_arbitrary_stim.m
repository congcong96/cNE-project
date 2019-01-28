function nedata = ne_calc_NE_with_arbitrary_stim(spk, trigger, binsize, plotopt)

if ~exist('plotopt','var')
    plotopt = 0;
end

% Get spike train matrix
spktrain = get_spktrain_from_arbitrary_trigger(spk, trigger, binsize);

if isfield(spk.spk, 'position')
    position = {spk.spk.position};
end

% Find cell assemblies
fprintf('\nDetecting Cell Assemblies\n');
[Patterns, Activities] = ca_detect_cell_assemblies_data(spktrain);

if plotopt == 1
    if exist('position','var')
        ca_plot_cell_assemblies_data(spktrain, Patterns, Activities, position);
        nedata.position = position;
    else
        ca_plot_cell_assemblies_data(spktrain, Patterns, Activities);
    end
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