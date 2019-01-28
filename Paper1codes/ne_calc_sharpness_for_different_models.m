function [pwc, sharpness] = ne_calc_sharpness_for_different_models(nefiles, varargin)

ip = inputParser;
addRequired(ip, 'nefiles', @iscell);
addParameter(ip, 'pwcinput', [], @(x) isempty(x) || isstruct(x));
addParameter(ip, 'numbins', 20, @(x) isscalar(x) && x > 0)
parse(ip, nefiles, varargin{:})

nefiles = ip.Results.nefiles;
pwcinput = ip.Results.pwcinput;
numbins = ip.Results.numbins;


for i = 1:length(nefiles)
    
    load(nefiles{i}, 'exp_site_nedata')
    basefile = regexp(nefiles{i}, '^\S+(?=(-ne))','match','once');

    NEneurons = ne_find_NE_pairs_or_groups(exp_site_nedata, 2);
    simspkfile = [basefile '-spktrainsforPWC.mat'];
    load(simspkfile);
    
    
    fsd = 200;
    stimdur = size(spktrain,2) * (1/fsd*1000);
    pval = 0.5;
    
    if isempty(pwcinput)
        pwc(i).real = pairwisecorr_function(spktrain,stimdur,fsd,NEneurons);
        [~, ~, sharpness(i).real, pd, ridx] = calc_pwc_sharpness_cdf(pwc.real, numbins, pval);
        pwc(i).real = pwc(i).real(ridx);
    else
        pwc = pwcinput;
        pd = [pwc(i).real.peakdelay];
        [~, ~, sharpness(i).real] = calc_pwc_sharpness_cdf(pwc.real, numbins, pval, 0, 1, pd);
    end
    
    
    
    if isempty(pwcinput)
        pwc(i).shuffled = pairwisecorr_function(shuffspktrain,stimdur,fsd,NEneurons);
        pwc(i).shuffled = pwc(i).shuffled(ridx);
    end
    [~, ~, sharpness(i).shuffled] = calc_pwc_sharpness_cdf(pwc.shuffled, numbins, pval, 0, 1, pd);   
    
    if isempty(pwcinput)
        pwc(i).PR = pairwisecorr_function(PRspktrain,stimdur,fsd,NEneurons);
        pwc(i).PR = pwc(i).PR(ridx);
    end
    [~, ~, sharpness(i).PR] = calc_pwc_sharpness_cdf(pwc.PR, numbins, pval, 0, 1, pd);
    
    
    if isempty(pwcinput)
        pwc(i).DG = pairwisecorr_function(DGspktrain,stimdur,fsd,NEneurons);
        pwc(i).DG = pwc(i).DG(ridx);
    end
    [~, ~, sharpness(i).DG] = calc_pwc_sharpness_cdf(pwc.DG, numbins, pval, 0, 1, pd);    


end

