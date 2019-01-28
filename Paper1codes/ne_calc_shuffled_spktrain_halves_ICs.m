function [ICshufflecell, varargout] = ne_calc_shuffled_spktrain_halves_ICs(spktrain,varargin)

% Get thresholds for cell assembly data to determine which units are
% 'considered part of an assembly'. Makes use of a random permutation test 
% to get the full range of ICA values by overriding the Marcenko-Pastur
% test and taking the first N PCs, where N = number of assemblies obtained
% from running the Marcenko-Pastur on the original dataset.
%
% Description of inputs:
%
%   spktrain:   Matrix of spike counts per time bin. Rows are single units;
%               columns are time bins.
%   perm:       If 'random', spike trains are shuffled randomly, preserving
%               only number of spikes in each spike train. If 'circular',
%               spike trains are circularly shifted, preserving the number
%               of spikes and all ISIs except for those at the edges.
%   num_iter:   Number of iterations for random permutation test. Default
%               set to 100.
%   prctile:    Confidence interval. For example, 95% returns 2.5 and 97.5
%               percentiles. Default set to [].
%   stdev:      Number of standard deviations from the mean. For example, 2
%               returns 2.5 and 97.5 percentile. Default set to 1.5. 
%   opts:       Struct array containing options for assembly_patterns.m.
%               Optional.
%
%   Default number of iterations for fast_ICA for original data is 500.
%   However, this is 100 for the random permutations for speed. This 
%   can be changed in opts. 
%
%   Updated 3/17/17 by js, added std for alpha

p = inputParser;
addOptional(p, 'perm', 'circular', @(x) strcmp(x, 'circular') || strcmp(x, 'random'))
addOptional(p, 'num_iter', 100, @(x) x >= 1)
addParameter(p, 'splitopt', 'interleaved', @(x) strcmp(x, 'interleaved') || strcmp(x,  'contiguous'))
addParameter(p, 'actopt', 0, @(x) x == 0 || x == 1)
addParameter(p, 'opts', [], @(x) isempty(x) || isstruct(x))
parse(p, varargin{:})

perm = p.Results.perm;
num_iter = p.Results.num_iter;
splitopt = p.Results.splitopt;
actopt = p.Results.actopt;
opts = p.Results.opts;

len = size(spktrain,2);

switch splitopt
    case 'interleaved'
        
        tenths = floor(len/10);
        idx = tenths:tenths:len;
        remainder = mod(len, 10);
        temp = zeros(1, 10);
        temp2 = 1:remainder;
        % split the remainder bins evenly
        temp(end-length(temp2)+1:end) = temp2;
        endidx = idx + temp;
        startidx = [1 endidx(1:end-1) + 1]; 
        segments = arrayfun(@(x,y) spktrain(:,x:y), startidx, endidx, 'UniformOutput',0);
        half1 = cell2mat(segments(1:2:end));
        half2 = cell2mat(segments(2:2:end));
        
    case 'contiguous'
        
        halflen = floor(len/2);
        half1 = spktrain(:,1:halflen);
        half2 = spktrain(:,halflen+1:end);
        
end

%process original data
if isempty(opts)
    patterns1 = assembly_patterns(half1);
    patterns2 = assembly_patterns(half2);
else
    patterns1 = assembly_patterns(half1, opts);
    patterns2 = assembly_patterns(half2, opts);
end

num_assemblies1 = size(patterns1,2);
num_assemblies2 = size(patterns2,2);

%for all the ICA values from permutation
ICshufflecell = cell(num_iter,2);
varargout{1} = cell(num_iter,1);

% zeros(size(patterns1,1), size(patterns1,2) + size(patterns2,2), num_iter);

for i = 1:num_iter
    
    fprintf('%d of %d iterations...\n',i,num_iter);
    spktrainperm = zeros(size(spktrain));
    [M, N] = size(spktrainperm);
    
    switch perm
        case 'random'
   
            for j = 1:M
                X = rand(N,1);
                [~,idx] = sort(X);
                spktrainperm(j,:) = spktrain(j,idx);      
            end

        case 'circular'
        
            shiftidx = random('unid',N,[M 1]);
    %         shiftidx = cell2mat(arrayfun(@(x) [x:N,1:x-1],shiftidx,'UniformOutput',0));            

            for j = 1:M
                spktrainperm(j,:) = circshift(spktrain(j,:),[0 shiftidx(j)]);
            end
    end
    
    switch splitopt
        case 'interleaved'
            segments = arrayfun(@(x,y) spktrainperm(:,x:y), startidx, endidx, 'UniformOutput',0);
            half1perm = cell2mat(segments(1:2:end));
            half2perm = cell2mat(segments(2:2:end));
        case 'contiguous'
            half1perm = spktrainperm(:,1:halflen);
            half2perm = spktrainperm(:,halflen+1:end);
    end
          
        
    %run each iteration through the modified assembly_patterns code
    if isempty(opts)
        ICshufflecell{i,1} = assembly_patterns_permutation(half1perm, ...
            num_assemblies1);
        ICshufflecell{i,2} = assembly_patterns_permutation(half2perm, ...
            num_assemblies2);
    else
        ICshufflecell{i,1} = assembly_patterns_permutation(half1perm, ...
            num_assemblies1, opts);
        ICshufflecell{i,2} = assembly_patterns_permutation(half2perm, ...
            num_assemblies2, opts);
    end
    
    if actopt == 1
        
        act1 = assembly_activity(ICshufflecell{i,1}, half2);
        act2 = assembly_activity(ICshufflecell{i,2}, half2);
        varargout{1}{i} = corr(act1',act2');
    end
        


    
end


    
