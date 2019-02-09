function shuffledICweights = ne_shuffle_IC_weights(exp_site_nedata, shuffleopt, varargin)

% shuffle IC weights to get control 'ensemble activity'

% initialize variables
nedata = exp_site_nedata.nedata;
ICweights = nedata.Patterns;
numneurons = size(ICweights,1);
NEmembers = nedata.NEmembers;
nonNEmembers = cellfun(@(x) setdiff(1:numneurons, x), NEmembers, 'UniformOutput',0);

if strcmp(shuffleopt, 'pseudorandom')
    newnonNEmembers = cellfun(@(x) setdiff(1:numneurons, x), varargin{1}, 'UniformOutput',0); 
end
% 
% if strcmp(shuffleopt, 'pseudorandom')
%     
%     assert(~isempty(varargin));
% %     NEmembers = nedata.NEmembers;
% %     NEsize = cellfun('length',NEmembers);
% %     NEsigneurons = cell(length(NEsize),1);
% %     
% %     for j = 1:length(NEsize)
% %         [NEsigneurons{j},~] = ne_find_non_NE_pairs_or_groups(exp_site_nedata, NEsize(j), 50, 0);
% %     end
% 
% end

shuffledICweights = zeros(size(ICweights));

for i = 1:size(ICweights,2)
    
    switch shuffleopt
        case 'random'
            shuffleidx = randsample(numneurons,numneurons);
            shuffledICweights(:,i) = ICweights(shuffleidx,i);
        case 'circular'
            if i == 1
                circidx = randsample(numneurons, size(ICweights,2));
            end
            shuffledICweights(:,i) = circshift(ICweights(:,i), circidx(i));
        case 'pseudorandom'
            sigweights = ICweights(NEmembers{i}, (i));
            sigidx = randperm(length(varargin{1}{i}));
            shuffledICweights(varargin{1}{i}(sigidx),i) = sigweights;
            nonsigweights = ICweights(nonNEmembers{i}, (i));
            nonsigidx = randperm(numneurons - length(varargin{1}{i}));
            shuffledICweights(newnonNEmembers{i}(nonsigidx),i) = nonsigweights;    
            
        otherwise
            error('Input #2 should either be ''random'', ''pseudorandom'' or ''circular''.')
    end

end