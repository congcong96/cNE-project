function fileids = ne_get_unique_file_ids(identifier)

if nargin == 0
    files = gfn('*dft.mat');
else
    files = gfn(sprintf('*%s*', identifier));
end
fileids = unique(cellfun(@(x) x(1:19), files, 'UniformOutput',0));


end

