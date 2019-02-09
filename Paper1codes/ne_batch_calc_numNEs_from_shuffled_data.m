function ne_batch_calc_numNEs_from_shuffled_data(files)

for i = 1:length(files)
    
    load(files{i})
    
    [~, ~, NEstats] = ne_calc_real_coincidence_with_models...
        (exp_site_nedata, 'sameNE', 1, 'preISI', 100, 1);
    
    id = regexp(files{i}, '^\S+(?=(-\d{3,4}um))','match','once');
    stim = regexp(files{i}, '(?<=(db-))rn\d{1,2}(?=(-fs))','match','once');
    dft = regexp(files{i}, '(?<=(-ne-))\d{1,3}dft(?=(.mat))','match','once');
    outfile = [id '-' stim '-' dft '-shuffledNEstats'];
    save(outfile,'NEstats')
    clear('NEstats')
end
    