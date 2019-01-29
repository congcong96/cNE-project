function ne_plot_constructive_MU_NE_STA_info_ptd(con_sta_MI_PTD, paramopt)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('paramopt', 'var')
    paramopt = 'nonparam';
end

NE_ptd = [con_sta_MI_PTD.NE_ptd];
MU_ptd = [con_sta_MI_PTD.MU_ptd];

NE_mi = [con_sta_MI_PTD.NE_info_extrap];
MU_mi = [con_sta_MI_PTD.MU_info_extrap];


% plot PTD graph
figure;
scatter(log10(MU_ptd), log10(NE_ptd), 10);
h = refline(1,0);
h.Color = 'r';
h.LineStyle = '--';
xlabel('Log_{10}(Multi-unit STRF PTD)')
ylabel('Log_{10}(cNE STRF PTD)')
tickpref;

switch paramopt
    case 'param'
        [~, p] = ttest(MU_ptd, NE_ptd);
    case 'nonparam'     
        p = signrank(MU_ptd, NE_ptd);
end
print_n_and_p(3/4, 1/4, length(MU_ptd), p); 
print_mfilename(mfilename);


% plot MI graph
figure;
scatter(log10(MU_mi), log10(NE_mi), 10);
h = refline(1, 0);
h.Color = 'r';
h.LineStyle = '--';
xlabel('Log_{10}(Multi-unit STRF MI)')
ylabel('Log_{10}(cNE STRF MI)')
tickpref;

switch paramopt
    case 'param'
        [~, p] = ttest(MU_mi, NE_mi);
    case 'nonparam'     
        p = signrank(MU_mi, NE_mi);
end
print_n_and_p(3/4, 1/4, length(MU_mi), p); 
print_mfilename(mfilename);

end

