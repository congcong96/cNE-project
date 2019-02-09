function ne_super_batch_save_upsampled_NEactivity(fileid, stim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(fileid)
    
    for j = 1:length(stim)
        
        [actstruct, usNE] = ne_batch_calc_upsampled_NEactivity(fileid{i}, stim{j});
        outfile = [fileid{i} '_' stim{j} '_upsampled_NE'];
        save(outfile, 'usNE', 'actstruct')
        clear('actstruct','usNE')
        
    end
end


end

