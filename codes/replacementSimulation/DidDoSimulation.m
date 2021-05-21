function res = DidDoSimulation(thistbl, roilabel, ntimes,  maxlen)

st = tic;
if nargin<3
    ntimes = 100;
    maxlen = 10;
end
if nargin<4
    maxlen = 10;
end
nsub = height(thistbl);
initmp = zeros(nsub, maxlen, ntimes);
res.after_corr = initmp;
res.after_auc = initmp;
res.after_aucbeh = initmp;
res.after_corruncer_raw = initmp;
res.after_corruncer_resid = initmp;
% res.after_corr_orth = initmp;

for jj = 1:nsub
    % read behavior
    clearidx = thistbl.resp{jj}~=-1 & thistbl.conf{jj}~=-1;
    diffc = thistbl.diff{jj}(clearidx);
    rtime = thistbl.rtime{jj}(clearidx);
    conf = thistbl.conf{jj}(clearidx);
    isco = thistbl.isc{jj}(clearidx);
    uncer = 5 - conf;
    
    % read Neural
    switch lower(roilabel)
        case 'dacc'
            roibeta_dmpfc_raw = mean(thistbl.DACC_JSH_oconf20_raw{jj},2);
            roibeta_rm = mean(thistbl.DACC_JSH_oconf20_rmDiffRT{jj},2);
        case 'dmpfc'
            roibeta_raw = mean(thistbl.DMPFC_JSH_oconf20_raw{jj},2);
            roibeta_rm = mean(thistbl.DMPFC_JSH_oconf20_rmDiffRT{jj},2);
        case 'fpc'
            roibeta_raw = mean(thistbl.FPCL26_JSH_oconf_raw{jj},2);
            roibeta_rm = mean(thistbl.FPCL26_JSH_oconf_rmDIFFRTONLY{jj},2);
    end
    
    if length(uncer) == length(roibeta_rm)
        ntrial = length(uncer);
    end
    % Before
    [~,~,stats] = glmfit([diffc, rtime], uncer);
    res.before_corr(jj,1) = corr(stats.resid, roibeta_rm);
    res.before_auc(jj,1) = fastAUC(isco*2-1, roibeta_rm, -1);
    
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
    
    res.before_aucbeh(jj,1) = fastAUC(sort_isco*2-1, sort_uncer_resid, -1);
    
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
            
            res.after_corr(jj, bb, tt) = corr(shuff_uncer_resid, sort_roibeta);
%             after_corruncer_raw(jj, bb, tt) = corr(shuff_uncer, sort_uncer(shufflefilter));
            res.after_corruncer_resid(jj, bb, tt) = corr(shuff_uncer_resid, sort_uncer_resid);
            res.after_auc(jj, bb, tt) = fastAUC(sort_isco*2-1, shuff_roibeta, -1);
            res.after_aucbeh(jj, bb, tt) = fastAUC(sort_isco*2-1, shuff_uncer_resid, -1);
%             res = Schmidt([sort_uncer_resid, shuff_uncer_resid]);
%             after_corr_orth(jj, bb, tt) = corr(res(:,2), sort_roibeta);
        end
    end
    
    
end


toc(st)

end