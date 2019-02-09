function varargout = ne_calc_sponevoked_intra_IC_correlations(exp_site_nedata)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nedata = exp_site_nedata.nedata;
c = 1;

if isfield(nedata, 'spon_patterns')
    ICspon = nedata.spon_patterns;
    sponscorrmat = corr(ICspon);
    varargout{c} = sponscorrmat(logical(triu(ones(size(ICspon,2)),1)));
    c = c+1;
end

if isfield(nedata,'evoked_patterns')
    ICevoked = nedata.evoked_patterns;
    evokedcorrmat = corr(ICevoked);
    varargout{c} = evokedcorrmat(logical(triu(ones(size(ICevoked,2)),1)));
    c = c+1;
end

if isfield(nedata,'Patterns')
    IC = nedata.Patterns;
    corrmat = corr(IC);
    varargout{c} = corrmat(logical(triu(ones(size(IC,2)),1)));
end


% sponevokedcorrmat = corr(ICspon, ICevoked);
% 
% maxsponcorr = max(sponevokedcorrmat, [], 2);
% maxevokedcorr = max(sponevokedcorrmat, [], 1);
% 
% 
% minval = min([sponcorr; evokedcorr; maxsponcorr; maxevokedcorr']);
% minval = floor(minval * 20) / 20;
% edges = minval:0.05:1;
% 
% figure;
% hold on;
% histogram(sponcorr, edges, 'Normalization', 'probability')
% histogram(evokedcorr, edges, 'Normalization', 'probability')
% histogram(maxsponcorr, edges, 'Normalization', 'probability')
% histogram(maxevokedcorr, edges, 'Normalization', 'probability')
% legend('Spontaneous', 'Evoked', 'Max spon-evoked', 'Max evoked-spon')


end

