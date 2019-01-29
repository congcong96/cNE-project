function NEsubset = ne_get_NEsubsets_from_member_spiketrain(NEtrain, membertrains, type)

% NEtrain: should be in STA resolution (i.e. 5 ms time bins)
% membertrains: can be in STA resolution or cNE resolution (i.e. 5 or 10 ms
% respectively)

if size(NEtrain, 2) == size(membertrains,2)
    membertrains = downsample_spiketrain(membertrains, 2);
end

if mod(size(NEtrain, 2), 2)
    NEtrain = [NEtrain 0];
    adjustopt = 1;
else
    adjustopt = 0;
end

NEmat = reshape(NEtrain, 2, []);
dsNEtrain = sum(NEmat, 1);

switch type
    case 'with/without'
        
        NEsubset = zeros(2, size(NEtrain, 2));

        comptrain = [dsNEtrain;sum(membertrains,1)];
        
        templogitrain = comptrain(1,:) == 1 & comptrain(2,:) >= 1;    
        tempNEtrain = zeros(size(NEmat));
        tempNEtrain(:, templogitrain) = NEmat(:,templogitrain);

        NEsubset(1,:) = reshape(tempNEtrain, 1, []);
        NEsubset(2,:) = NEtrain - NEsubset(1,:);
        
    case 'individual_members'
        
        NEsubset = zeros(size(membertrains,1), size(NEtrain, 2));
        
        for i = 1:size(membertrains, 1)
    
            comptrain = [dsNEtrain;membertrains(i,:)];
            
            templogitrain = comptrain(1,:) == 1 & comptrain(2,:) >= 1;
            tempNEtrain = zeros(size(NEmat));
            tempNEtrain(:, templogitrain) = NEmat(:,templogitrain);

            NEsubset(i,:) = reshape(tempNEtrain, 1, []);

        end        
end

if adjustopt
    NEsubset = NEsubset(:,1:end-1);
end

