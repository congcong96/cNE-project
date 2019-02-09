function actstruct = ne_project_upsampled_activities_for_diff_binsizes(ICweights, spktrain)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fn = fieldnames(ICweights);

for i = 1:length(fn)
    
    fprintf('\nProcessing %d of %d bin sizes...', i, length(fn))
    
    patterns = ICweights.(fn{i});
    
    for j = 1:size(patterns,2)
        actstruct.(fn{i})(j,:) = assembly_activity(patterns(:,j), spktrain);
    end
    

end

fprintf('\n')


