function ne_batch_calc_spon_cNEs(files, DF)

for i = 1:length(files)
    
    load(files{i})
    
    if exist('spon_spktrain', 'var')
        spktrain = spon_spktrain;
        clear('spon_spktrain','spon_edges');
    end
    
    ds_spktrain = downsample_spiketrain(spktrain, DF);
    
    % Find cell assemblies
    fprintf('\nDetecting Cell Assemblies\n');
    [Patterns, Activities] = ca_detect_cell_assemblies_data(ds_spktrain);
    
    % Assign data to struct for output argument
    nedata.spktrain = ds_spktrain;
    nedata.fsdvd = 96000;
    nedata.df = DF;
    nedata.position = {spk.spk.position};
    nedata.Patterns = Patterns;
    nedata.Activities = Activities;
    
    parts = regexp(files{i}, '(?<part1>^\S+)-spon\w+-(?<part2>\S+).mat', 'names');
    outfile = [parts.part1 '-spon-' parts.part2 sprintf('-ne-%ddft.mat', DF)];
    save(outfile, 'nedata');
    
    exp_site_nedata = ne_create_exp_site_nedata_file(outfile);
    save(outfile, 'exp_site_nedata')
    

end

