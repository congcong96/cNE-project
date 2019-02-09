function [actstruct, usNE] = ne_batch_calc_upsampled_NEactivity(fileid, stim, dfs, reqdf)
% (spkfile, df, reqdf, varargin)
if ~exist('reqdf', 'var')
    reqdf = 2;
end

if ~exist('saveopt', 'var')
    saveopt = 0;
end

if ~exist('dfs','var') || isempty(dfs)
    dfs = [4 10 16 20 30 40 70 100];
end

% get stimfile
% for i = 1:length(stimulus)
    
% stim = stimulus{i};

drive = gcdr;
subfolder = 'Ripple_Noise\downsampled';
stimfolder = fullfile(drive,subfolder);

stimfilename = sprintf('%s-*_DFt1_DFf5_matrix.mat',stim);
stimfile = fullfile(stimfolder,stimfilename);
stimfile = dir(stimfile);
load(stimfile.name);

subfolder = 'Cell_Assemblies\PaperNEs';
spkfilepath = fullfile(drive, subfolder, fileid);
spkfile = gfn([spkfilepath '*' stim '-*' 'strfcmb.mat']);

assert(length(spkfile) == 1)

for j = 1:length(dfs)

    usNE{j} = ne_calc_offset_NEactivity(spkfile{1}, dfs(j), reqdf, stim_mat);
    temp = zeros(size(usNE{j}(1).act,1), size(usNE{j}(1).act,2) * length(usNE{j}));

    for k = 1:size(usNE{j}(1).act,1)            
        for ii = 1:length(usNE{j})                
            temp(k, ii:length(usNE{j}):end) = usNE{j}(ii).act(k,:);
        end
    end

    actstruct.(sprintf('df%d',dfs(j))) = temp;
end

% if saveopt == 1
%     outfile = [fileid '_' stim '_upsampled_NE'];
%     save(outfile, 'usNE', 'actstruct')
% end

% end
end

