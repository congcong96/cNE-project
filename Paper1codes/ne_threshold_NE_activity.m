function [actbin, threshold] = ne_threshold_NE_activity(actstruct, threshopt, varargin)

if nargin == 1
    threshopt = 'mad';
end

fn = fieldnames(actstruct);

switch threshopt
    case 'fixednum'
        
        for i = 1:length(fn)

            actmat = actstruct.(fn{i});
            tempactbin = zeros(size(actmat));

            [b, idx] = sort(actmat,2, 'descend');

            for j = 1:size(idx,1)
                tempactbin(j,idx(j,1:varargin{1})) = 1;
            end

            threshold.(fn{i}) = b(:,varargin{1});

            actbin.(fn{i}) = tempactbin;

        end
        
    case 'mad'
        
        for i = 1:length(fn)
            
            actmat = actstruct.(fn{i});
            
            dev = mad(actmat,0,2);
            m = mean(actmat,2);
            threshold.(fn{i}) = m + 5*dev;
            
            for j = 1:size(actmat,1)
                actbin.(fn{i})(j,:) = actmat(j,:) >= threshold.(fn{i})(j);
            end
        end
        
end