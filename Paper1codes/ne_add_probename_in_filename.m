function ne_add_probename_in_filename(files, probename)

% to be used for spkfiles to add probename in filename for easy processing
% Written 7/6/17 by JS

for i = 1:length(files)
    
    parts = regexp(files{i}, '(?<=(\d{1,2}))-(?=(fs))', 'split');
    assert(length(parts) == 2);
    newname = [parts{1} '-' probename '-' parts{2}];
    movefile(files{i}, newname);

end

