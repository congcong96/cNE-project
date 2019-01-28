function [corrmat, IC1, IC2] = ne_compare_ICweights_of_halves(exp_site_nedata, splitopt)

if ~exist('splitopt','var')
    splitopt = 'interleaved';
end

dsspktrain = exp_site_nedata.nedata.spktrain;

switch splitopt
    case 'interleaved'
        % split the remainder bins evenly
        len = size(dsspktrain,2);
        tenths = floor(len/10);
        idx = tenths:tenths:len;
        remainder = mod(len, 10);
        temp = zeros(1, 10);
        temp2 = 1:remainder;
        temp(end-length(temp2)+1:end) = temp2;
        endidx = idx + temp;
        startidx = [1 endidx(1:end-1) + 1]; 
        segments = arrayfun(@(x,y) dsspktrain(:,x:y), startidx, endidx, ...
            'UniformOutput',0);
        half1 = cell2mat(segments(1:2:end));
        half2 = cell2mat(segments(2:2:end));

%         for i = 1:length(idx)
%             if mod(i,2) %odd
% 
%                 if i == 1
%                     half1 = [half1 dsspktrain(:,1:idx(i))];
%                 else
%                     half1 = [half1 dsspktrain(:, idx(i-1) + 1:idx(i))];
%                 end
% 
%             else %even
% 
%                 half2 = [half2 dsspktrain(:, idx(i-1) + 1 : idx(i))];
% 
%             end
%         end
        
    case 'contiguous'
        len = size(dsspktrain,2);
        halves = floor(len/2);
        half1 = dsspktrain(:,1:halves);
        half2 = dsspktrain(:,halves+1:end);
end
        
IC1 = assembly_patterns(half1);
IC2 = assembly_patterns(half2);

close all

corrmat = abs(corr(IC1, IC2));


end

