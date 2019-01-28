function ne_add_stimlength_in_filename(files, stimlength)

% to be used for spkfiles to add stimlength in filename for easy processing
% Written 7/6/17 by JS

for i = 1:length(files)
    
    parts = regexp(files{i}, '(?<=(\d{1,2}))-(?=(a1))', 'split');
    assert(length(parts) == 2);
    if ischar(stimlength)
        newname = [parts{1} '-' stimlength '-' parts{2}];
    else
        assert(length(stimlength) == length(files))
        newname = [parts{1} '-' stimlength(i) '-' parts{2}];
    end
    movefile(files{i}, newname);

end
