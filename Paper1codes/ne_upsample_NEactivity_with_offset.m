function [NEraster, usNE, temp, NEthresh] = ne_upsample_NEactivity_with_offset(spkfile, df, reqdf, stimmat)

usNE = ne_calc_offset_NEactivity(spkfile, df, reqdf, 1, stimmat);    
NEthresh = mean(cell2mat({usNE.NEthresh}), 2);

temp = zeros(size(usNE(1).Activities,1), size(usNE(1).Activities,2) * 2);

for i = 1:length(usNE)        
    temp(:,i:length(usNE):end) = usNE(i).Activities;        
end

NEraster = zeros(size(temp));
for j = 1:size(temp,1)        
    NEraster(j,:) = temp(j,:) >= NEthresh(j);
end