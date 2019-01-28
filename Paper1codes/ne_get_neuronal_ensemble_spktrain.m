function NE_train = ne_get_neuronal_ensemble_spktrain(exp_site_nedata, varargin)

% Gets cNE logical with member neuron spktrains/logicals.
%
%   exp_site_nedata: standard input from initial cell assembly analysis.
%   
%   memeuopt: Member neuron option. Include member neurons below cNE train
%   if desired. 'none' means that they are not returned; 'logical' returns
%   binary neuronal spike trains; and 'raw' returns neuronal spike trains as
%   they are.
%
%   threshalpha: Threshold chosen for cNE activity. Takes values from 99.0
%   to 99.9 in 0.1 increments.
%   
%   method: Takes 'repeat' or 'interpolate'. 'repeat' first thresholds
%   cNE activity and upsamples by repeating value for bin groups. Best used
%   for splitting neuronal spike groups. 'Interpolate' upsamples by
%   interpolating first before thresholding. 
%
%   Jermyn See, updated 4/6/18. Cleaned up code and made it more efficient 
%   and logical.


ip = inputParser;
addRequired(ip, 'exp_site_nedata', @isstruct)
addParameter(ip, 'memneuopt', 'raw', @(x) ischar(x) && any(ismember(x, {'none','logical','raw'})))
addParameter(ip, 'threshalpha', 99.5, @(x) x >= 99 && x <= 99.9)
addParameter(ip, 'method', 'repeat', @(x) strcmp(x, 'repeat') || strcmp(x, 'interpolate'))
parse(ip, exp_site_nedata, varargin{:})

exp_site_nedata = ip.Results.exp_site_nedata;
memneuopt = ip.Results.memneuopt;
threshalpha = ip.Results.threshalpha;
method = ip.Results.method;

nedata = exp_site_nedata.nedata;
stadf = 10;
df = exp_site_nedata.df;
assert(stadf <= df);

NEthresh_alpha = nedata.NEthresh_alpha;
NEthresh = nedata.NEthresh;
NEact = nedata.Activities;

alphaidx = NEthresh_alpha == threshalpha;
NEthresh = NEthresh(alphaidx, :);

usfactor = df/stadf;
NE_train = zeros(size(NEact,1), size(NEact,2) * usfactor);

switch method
    
    case 'repeat'
        
        for i = 1:length(NEthresh)            
            temptrain = NEact(i,:) >= NEthresh(i);            
            for j = 1:usfactor
                NE_train(i,j:usfactor:end) = temptrain;
            end
        end
        
    case 'interpolate'
        
        for i = 1:length(NEthresh)            
            temptrain = interp(NEact, usfactor);
            NE_train(:,i) = temptrain >= NEthresh(i);
        end
end

sta_spktrain = nedata.sta_spktrain;
sizediff = size(NE_train,2) - size(sta_spktrain,2);
NE_train = NE_train(:, 1:end-sizediff);

if ~strcmp(memneuopt, 'none')
    
    NEmembers = nedata.NEmembers;
    
    if strcmp(memneuopt, 'logical')       
        sta_spktrain = sta_spktrain >= 1;        
    end
        
    NE_train = cellfun(@(x,y) [x; sta_spktrain(y,:)], mat2cell...
        (NE_train, ones(size(NE_train,1),1), size(NE_train,2)),...
        NEmembers, 'UniformOutput', 0);
        
end    




