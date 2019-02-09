function [corrvals, varargout] = ne_calc_sta_corrcoef_between_real_and_surr_spiketrains(exp_site_nedata, normopt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nargoutchk(1,2)

if normopt == 1
    stamat = ne_calc_sig_sta(exp_site_nedata);

    for i = 1:4

        [~,staPred(i).preRF,corrvals(i).preRF] = ne_get_surrspktrain_projval_bin_shuffle(exp_site_nedata,...
            'random', 15, stamat, 1);
        [~,staPred(i).prePFRRF,corrvals(i).prePFRRF] = ne_get_surrspktrain_projval_bin_shuffle(exp_site_nedata,...
            'frprob', 15, stamat, 1);

    end
    
else
    
    for i = 1:1
        
        spktrainPred = ne_circularly_shuffle_spkmatrix (exp_site_nedata, 'sta');
        stim = regexp(exp_site_nedata.stim, '^rn\d{1,2}','match','once');
        load(sprintf('I:\\Ripple_Noise\\downsampled_for_MID\\%s-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt1_DFf5_matrix.mat', stim))
        stim_mat = stim_mat(:,1:10:end);
        staPred(i).shuffled = quick_calc_sta(stim_mat, spktrainPred, 20, 64);
        corrvals(i).shuffled = diag(corr(staPred(i).shuffled', exp_site_nedata.nedata.stamat'));
        
        spktrainPred = ne_generate_dg_simulated_spktrains(exp_site_nedata.nedata.sta_spktrain);
        staPred(i).DG = quick_calc_sta(stim_mat, spktrainPred, 20, 64);
        corrvals(i).DG = diag(corr(staPred(i).DG', exp_site_nedata.nedata.stamat'));

        
        [~,staPred(i).PR, corrvals(i).PR] = ne_get_surrspktrain_projval_bin_shuffle(exp_site_nedata,...
            'random');
%         [~,staPred(i).prePFRRF,corrvals(i).prePFRRF] = ne_get_surrspktrain_projval_bin_shuffle(exp_site_nedata,...
%             'frprob');
        

    end
    
end
    
varargout{1} = staPred; 
    
return



