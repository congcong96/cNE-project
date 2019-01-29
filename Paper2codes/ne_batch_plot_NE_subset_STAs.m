function ne_batch_plot_NE_subset_STAs(sig_sta, sigopt, type, saveopt)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('sigopt', 'var')
    sigopt = 0;
end

if ~exist('type','var')
    type = 'with/without';
end

if ~exist('saveopt', 'var')
    saveopt = [];
end

assert(isempty(saveopt) || ischar(saveopt))

if ~isempty(saveopt)
    close all
    pdffiles = cell(length(sig_sta), 1);
end

for i = 1:length(sig_sta)
    
    load(sig_sta(i).filename, 'exp_site_nedata');
    ne_plot_NE_subset_STAs(exp_site_nedata, sig_sta(i).NE, sigopt, type);
%     suptitle(sprintf('%s_NE%d', sig_sta(i).filename, sig_sta(i).NE));
    
    if ~isempty(saveopt)
        pdffiles{i} = sprintf('%s%d.pdf', saveopt, i);
        export_fig(gcf, pdffiles{i}, '-nocrop'); 
        close all
    end

end

if saveopt
    append_pdfs(sprintf('%s.pdf', saveopt), pdffiles{:});
    cellfun(@delete, pdffiles);
end

