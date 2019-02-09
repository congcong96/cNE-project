function ne_batch_calc_NEsta_info_for_binsizes(fileids, df, stimtype, USopt)


if nargin == 2
    stimtype = {'rn1','rn4','rn8','rn16'};
end

if ischar(stimtype)
    stimtype = {stimtype};
end

for i = 1:length(fileids)
    
    sitenum = regexp(fileids{i},'site\d{1,2}$','match','once');
    date = regexp(fileids{i}, '^\d{6}(?=(_\d{6}))','match','once');
    
    for j = 1:length(stimtype)        
        
        fprintf('\nProcessing %s_%s\n', sitenum, stimtype{j})

        info = ne_calc_NEsta_info_for_binsizes(fileids{i}, stimtype{j}, df, USopt);

        clc;
        
        save([sprintf('%s_%s_%s',date,sitenum,stimtype{j}) '_info_offset'], 'info')
        clear('info');
         
    end
    
end

% save('NEstainfo_vs_binsize_interpolate2','infoproj')

return
        
