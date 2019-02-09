function [npat1, npat2] = ne_arrange_ICweight_corrmat(pat1, pat2, thresh, varargin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
    CI1 = varargin{1};
    CI2 = varargin{2};
    pat1 = pat1 <= CI1(1) | pat1 >= CI1(2);
    pat2 = pat2 <= CI2(1) | pat2 >= CI2(2);
end


corrmat = abs(corr(pat1,pat2));

dim1 = size(corrmat,1);
dim2 = size(corrmat,2);
mindim = min([dim1 dim2]);

sigcorrvals = sort(corrmat(corrmat>thresh), 'descend');

npat1 = zeros(size(pat1));
npat2 = zeros(size(pat2));

list1 = 1:dim1;
list2 = 1:dim2;

c = 1;

for i = 1:length(sigcorrvals)
    
    [x,y] = find(sigcorrvals(i) == corrmat);
    
    if ismember(x,list1) && ismember(y,list2) && c <= mindim
        list1 = setdiff(list1, x);
        list2 = setdiff(list2, y);
           
        npat1(:,c) = pat1(:,x);
        npat2(:,c) = pat2(:,y);
        
        c = c+1;        
    else
        continue
    end


end

c1 = c;

for i = 1:length(list1)    
    npat1(:, c1) = pat1(:, list1(i));
    c1 = c1+1;
end

c2 = c;

for i = 1:length(list2)
    npat2(:, c2) = pat2(:, list2(i));
    c2 = c2+1;
end



