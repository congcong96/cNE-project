function dist = ne_calc_interneuronal_distances(comb, pos, distopt)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

dist = zeros(size(comb,1),1);

switch distopt
    case 'channel'
        horidist = min(diff(unique(pos(:,1))));
        vertdist = min(diff(unique(pos(:,2))));
end

for i = 1:size(comb,1)
    pos1 = pos(comb(i,1),:);
    pos2 = pos(comb(i,2),:);
    temp = pos1 - pos2;
    
    switch distopt
        case 'distance'    
            dist(i) = hypot(temp(1), temp(2));
        case 'channel'
            dist(i) = max(abs([temp(1)./horidist temp(2)./vertdist]));
    end
end


end

