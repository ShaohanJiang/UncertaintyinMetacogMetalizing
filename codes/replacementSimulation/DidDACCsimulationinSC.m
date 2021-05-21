%% DACC in SC
clear; home;
load('../BehavNeuralTbl.mat', 'BehavNeuraltbl_SC')
st = tic;
thistbl = BehavNeuraltbl_SC;
nsub = height(thistbl);
maxlen = 10;
ntimes = 1e4;
initmp = zeros(nsub, maxlen, ntimes);
after_corr = initmp;
after_auc = initmp;
after_aucbeh = initmp;
after_corruncer_raw = initmp;
after_corruncer_resid = initmp;
% after_corr_orth = initmp;

for jj = 1:nsub
    % behavior
    clearidx = thistbl.resp{jj}~=-1 & thistbl.conf{jj}~=-1;
    diffc = thistbl.diff{jj}(clearidx);
    rtime = thistbl.rtime{jj}(clearidx);
    conf = thistbl.conf{jj}(clearidx);
    isco = thistbl.isc{jj}(clearidx);
    uncer = 5 - conf;
    % Neural
    roibeta_dmpfc_raw = mean(thistbl.DACC_JSH_oconf20_raw{jj},2);
    roibeta_rm = mean(thistbl.DACC_JSH_oconf20_rmDiffRT{jj},2);
    
    if length(uncer) == length(roibeta_rm)
        ntrial = length(uncer);
    end
    % Before
    [~,~,stats] = glmfit([diffc, rtime], uncer);
    before_corr(jj,1) = corr(stats.resid, roibeta_rm);
    before_auc(jj,1) = fastAUC(isco*2-1, roibeta_rm, -1);
    
    tmpmatrix = [uncer, diffc, rtime, (diffc-mean(diffc)).*(rtime-median(rtime)), roibeta_rm, isco];
    tmpmatrix = sortrows(tmpmatrix, 4);
    sort_roibeta = tmpmatrix(:, end-1);
    sort_isco = tmpmatrix(:, end);
    sort_uncer = tmpmatrix(:, 1);
    sort_diffrt = tmpmatrix(:, 2:3);
    
    dmat = [ones(length(diffc),1), sort_diffrt];
    idmat = pinv(dmat);
    beta = idmat*sort_uncer;
    sort_uncer_resid = sort_uncer - dmat*beta;
    
    before_aucbeh(jj,1) = fastAUC(sort_isco*2-1, sort_uncer_resid, -1);
    
    for bb = 1:maxlen
        tmpshuffle = zeros(size(tmpmatrix));
        % shuffle
        nbin = ceil(ntrial/bb);
%         binlabel = repmat(1:nbin, bb, 1);
%         binlabel = binlabel(1:ntrial)';
        for tt = 1:ntimes
            
            for ii = 1:nbin
                if ii~=nbin
                    thisidx = bb*(ii-1)+1 : bb*ii;
                else
                    thisidx = bb*(ii-1)+1 : ntrial;
                end
                thisbin = tmpmatrix(thisidx, :);
                for mm = 1:length(thisidx)
                    tmpshuffle(thisidx(mm), :) = thisbin(randperm(length(thisidx),1), :);
                end
            end

            shuff_uncer = tmpshuffle(:, 1);
            shuff_diffrt = tmpshuffle(:, 2:3);
            shuff_roibeta = tmpshuffle(:, end-1);
            
            dmat = [ones(length(shuff_diffrt),1), shuff_diffrt];
            idmat = pinv(dmat);
            beta = idmat*shuff_uncer;
            shuff_uncer_resid = shuff_uncer - dmat*beta;
            
            after_corr(jj, bb, tt) = corr(shuff_uncer_resid, sort_roibeta);
%             after_corruncer_raw(jj, bb, tt) = corr(shuff_uncer, sort_uncer(shufflefilter));
            after_corruncer_resid(jj, bb, tt) = corr(shuff_uncer_resid, sort_uncer_resid);
            after_auc(jj, bb, tt) = fastAUC(sort_isco*2-1, shuff_roibeta, -1);
            after_aucbeh(jj, bb, tt) = fastAUC(sort_isco*2-1, shuff_uncer_resid, -1);
%             res = Schmidt([sort_uncer_resid, shuff_uncer_resid]);
%             after_corr_orth(jj, bb, tt) = corr(res(:,2), sort_roibeta);
        end
    end
    
    
end


toc(st)
% SAVE
save res_daccinsc_simulation_replacement_10k after* before* maxlen


%%

after_corrauc = zeros(ntimes, maxlen);

    for jj = 1:ntimes
        after_corrauc(jj, :) = corr(before_auc, after_auc(:,:,jj));

    end

figure; histogram(after_corrauc(:,10))
box off
hold on; plot([0,0], [0,700], '--', 'LineWidth',1.5)

