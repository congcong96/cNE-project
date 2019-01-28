function [spktrain, position, binedges] = ne_create_spktrain_matrix_wo_stim(spk, binsize, numbins)

if nargin == 2
    numbins = [];
end

%get max spktime

if ~isfield(spk, 'spiketimes')
    spk = spk.spk;
end

if isempty(numbins)
    
    try
%          allspk = cell2mat({spk.filt_spiketimes});
        allspk = cell2mat({spk.spiketimes});
    catch
        allspk = cell2mat({spk.spiketimes}');
    end
    maxspk = ceil(max(allspk));
    remainder = mod(ceil(maxspk),binsize);
    lastbinedge = (maxspk+remainder);
    
else
    lastbinedge = binsize * numbins;
end


binedges = 0:binsize:lastbinedge;

spktrain = zeros(size(spk,2),length(binedges)-1);
position = cell(length(spk),1);

for i = 1:length(spk)
    
spktimes = spk(i).spiketimes;
% keyboard
try
    spktrain(i,:) = histcounts(spktimes,binedges);
catch
    spktrain(i,:) = histc(spktimes,binedges(1:end-1));
end

if isfield(spk, 'position')
    position{i} = (spk(i).position);
end


end

