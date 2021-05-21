%%
function [res, getglm] = calcBasicDACC(thistbl)
st = tic;
% thistbl = BehavNeuraltbl_SC;
nsub = height(thistbl);
% figure;
for jj = 1:nsub
    % behavior
    clearidx = thistbl.resp{jj}~=-1 & thistbl.conf{jj}~=-1;
    diffc = thistbl.diff{jj}(clearidx);
    rtime = thistbl.rtime{jj}(clearidx);
    ctime = thistbl.ctime{jj}(clearidx);
    conf = thistbl.conf{jj}(clearidx);
    isco = thistbl.isc{jj}(clearidx);
    uncer = 5 - conf;
    absuncer = -1*abs(uncer-mean(uncer));
    
    getglm(jj).subnames = thistbl.subnames{jj};
    atmp = fitglm(zscore([diffc, rtime]), zscore(uncer));
    getglm(jj).betadiffc = atmp.Coefficients.Estimate(2);
    getglm(jj).betartime = atmp.Coefficients.Estimate(3);
    getglm(jj).residuncer = atmp.Residuals.Raw;
    getglm(jj).fituncer = atmp.Fitted.Response;
    getglm(jj).zscorevariables = atmp.Variables;
    getglm(jj).rawvariables = [diffc, rtime, uncer];
    
   
    
    dmat = [ones(length(diffc),1), zscore([diffc, rtime])];
    idmat = pinv(dmat);
    beta = idmat * zscore(uncer);
    fituncer = dmat*beta;
    uncer_resid = zscore(uncer) - fituncer;
    
    B = glmfit(zscore([diffc, rtime]), zscore(uncer));
    res.betadiffc(jj,1) = B(2);
    res.betartime(jj,1) = B(3);
    
    res.corr_absresidfituncer(jj,:) = [corr(abs(uncer_resid), fituncer), corr(abs(uncer_resid), uncer_resid)]; 
    res.corr_ctrt(jj) = corr(ctime, rtime);
    % ctime shape in different uncer rating
    ctime_bin{1} = ctime(uncer==1);
    ctime_bin{2} = ctime(uncer==2);
    ctime_bin{3} = ctime(uncer==3);
    ctime_bin{4} = ctime(uncer==4);
    res.ctime_bin(jj,:) = [mean(ctime_bin{1}),mean(ctime_bin{2}),mean(ctime_bin{3}),mean(ctime_bin{4})];
    
    bin_seed = rtime;
    bin_prctile =  prctile(bin_seed,[25,50,75]);
    nbin = length(bin_prctile)+1;
    
    for ii = 1:nbin
        if ii == 1
            res.uncer_resid_bin{jj,1} = uncer_resid(bin_seed<=bin_prctile(1));
            res.ctime_binbyrtime{jj,1} = ctime(bin_seed<=bin_prctile(1));
        elseif ii == nbin
            res.uncer_resid_bin{jj, nbin} = uncer_resid(bin_seed>bin_prctile(nbin-1));
            res.ctime_binbyrtime{jj, nbin} = ctime(bin_seed>bin_prctile(nbin-1));
        else
            res.uncer_resid_bin{jj,ii} = uncer_resid(bin_seed>bin_prctile(ii-1) & bin_seed<=bin_prctile(ii));
            res.ctime_binbyrtime{jj,ii} = ctime(bin_seed>bin_prctile(ii-1) & bin_seed<=bin_prctile(ii));
        end
        
%         res.uncer_resid_bin{jj,2} = uncer_resid(bin_seed>bin_prctile(1) & bin_seed<=bin_prctile(2));
%         res.uncer_resid_bin{jj,3} = uncer_resid(bin_seed>bin_prctile(2) & bin_seed<=bin_prctile(3));
%         res.uncer_resid_bin{jj,4} = uncer_resid(bin_seed>bin_prctile(3) & bin_seed<=bin_prctile(4));
%         res.uncer_resid_bin{jj,5} = uncer_resid(bin_seed>bin_prctile(4));
    end
    
    % Neural
    daccbeta_raw = mean(thistbl.DACC_JSH_oconf20_raw{jj},2); 
    daccbeta_rm = mean(thistbl.DACC_JSH_oconf20_rmMEANDIFFRT{jj},2); 
    
    dmpfcbeta_raw = mean(thistbl.DMPFC_JSH_oconf20_raw{jj},2);  
    dmpfcbeta_rm = mean(thistbl.DMPFC_JSH_oconf20_rmMEANDIFFRT{jj},2);
    
    fpcbeta_raw = mean(thistbl.FPCL26_JSH_oconf_raw{jj},2);  % fpc_l_6mm_oconf20_raw FPCL26_JSH_oconf_raw
    fpcbeta_rm = mean(thistbl.FPC_JSH_oconf20_rmMEANDIFFRT{jj},2); % FPCL26_JSH_oconf_rmDIFFRTONLY FPC_JSH_oconf20_rmMEANDIFFRT
    
    iplbeta_raw =  mean(thistbl.IPL_JSH_oconf20_raw{jj},2); %ipl_r_6mm_oconf20_raw IPL_JSH_oconf20_raw IPL_overlap_JSH_oconf20_raw
    iplbeta_rm = mean(thistbl.IPL_JSH_oconf20_rmDiffRT{jj},2);
    
    ifjbeta_raw = mean(thistbl.IFJ_conj_JSH_oconf20_raw{jj},2);  
    ifjbeta_rm = mean(thistbl.IFJ_conj_JSH_oconf20_rmMEANDIFFRT{jj},2);
    
    tpjbeta_raw = mean(thistbl.TPJ_JSH_oconf20_raw{jj},2);  
    tpjbeta_rm = mean(thistbl.TPJ_JSH_oconf20_rmMEANDIFFRT{jj},2);
    %-------------------------------------------------------------%
%     res.corr_ipldiffc_raw(jj,1) = corr(iplbeta_raw, diffc);
%     res.corr_iplrtime_raw(jj,1) = corr(iplbeta_raw, rtime);
%     res.corr_iplresiduncer_raw(jj,1) = corr(iplbeta_raw, uncer_resid);

%     res.corr_ipldiffc_rm(jj,1) = corr(iplbeta_rm, diffc);
%     res.corr_iplrtime_rm(jj,1) = corr(iplbeta_rm, rtime);
%     res.corr_iplresiduncer_rm(jj,1) = corr(iplbeta_rm, uncer_resid);
%     
%     res.corr_dmpfcdiffc_raw(jj,1) = corr(dmpfcbeta_raw, diffc);
%     res.corr_dmpfcrtime_raw(jj,1) = corr(dmpfcbeta_raw, rtime);
%     res.corr_dmpfcresiduncer_raw(jj,1) = corr(dmpfcbeta_raw, uncer_resid);
%     
%     res.corr_dmpfcdiffc_rm(jj,1) = corr(dmpfcbeta_rm, diffc);
%     res.corr_dmpfcrtime_rm(jj,1) = corr(dmpfcbeta_rm, rtime);
%     res.corr_dmpfcresiduncer_rm(jj,1) = corr(dmpfcbeta_rm, uncer_resid);
%     
%     res.corr_fpcdiffc_raw(jj,1) = corr(fpcbeta_raw, diffc);
%     res.corr_fpcrtime_raw(jj,1) = corr(fpcbeta_raw, rtime);
%     res.corr_fpcresiduncer_raw(jj,1) = corr(fpcbeta_raw, uncer_resid);
%     
%     res.corr_fpcdiffc_rm(jj,1) = corr(fpcbeta_rm, diffc);
%     res.corr_fpcrtime_rm(jj,1) = corr(fpcbeta_rm, rtime);
%     res.corr_fpcresiduncer_rm(jj,1) = corr(fpcbeta_rm, uncer_resid);
    
    orth_matrix = Schmidt([diffc, rtime, uncer_resid]);
    
    res.beta_ipl_diffcrtresiduncer_raw(jj,:) = glmfit(zscore(orth_matrix), zscore(iplbeta_raw), []);
    res.beta_dacc_diffcrtresiduncer_raw(jj,:) = glmfit(zscore(orth_matrix), zscore(daccbeta_raw), []);
    res.beta_dmpfc_diffcrtresiduncer_raw(jj,:) = glmfit(zscore(orth_matrix), zscore(dmpfcbeta_raw), []);
    res.beta_fpc_diffcrtresiduncer_raw(jj,:) = glmfit(zscore(orth_matrix), zscore(fpcbeta_raw), []);
    
    
    %-------------------------------------------------------------------%
    
    roibeta_raw = daccbeta_raw;
    roibeta_rm = daccbeta_rm;
    
    if length(uncer) == length(daccbeta_rm)
        ntrial = length(uncer);
    else
        error('Error incur!')
    end
    
    
    % Before
    [~,~,stats] = glmfit([diffc, rtime], uncer);
    res.before_corr(jj,1) = corr(stats.resid, roibeta_rm);
    res.before_corr_raw(jj,1) = corr(stats.resid, roibeta_raw);
    res.before_aucbeh1(jj,1) = fastAUC(isco*2-1, stats.resid, -1);
    res.before_auc(jj,1) = fastAUC(isco*2-1, roibeta_rm, -1);
    res.before_aucbeh2(jj,1) = fastAUC(isco*2-1, uncer_resid, -1);
    
    res.roi_auc(jj,:) = [fastAUC(isco*2-1, daccbeta_rm, -1), fastAUC(isco*2-1, dmpfcbeta_rm, -1),...
        fastAUC(isco*2-1, fpcbeta_rm, -1), fastAUC(isco*2-1, ifjbeta_rm, -1), fastAUC(isco*2-1, tpjbeta_rm, -1)];
    
    tmpmatrix = [uncer, diffc, rtime, (diffc-mean(diffc)).*(rtime-median(rtime)), ctime, fituncer, ...
        uncer_resid, roibeta_rm, isco];
    tmpmatrix = sortrows(tmpmatrix, 6);
    
    sort_neghf = tmpmatrix(tmpmatrix(:,7)<0,:);
    sort_poshf = tmpmatrix(tmpmatrix(:,7)>0,:);
    
    sort_roibeta = tmpmatrix(:, end-1);
    sort_isco = tmpmatrix(:, end);
    sort_uncer = tmpmatrix(:, 1);
    sort_diffrt = tmpmatrix(:, 2:3);
    sort_uncer_resid = tmpmatrix(:, end-2);
    sort_fitresid = tmpmatrix(:, end-3);
    
    res.corr_roibetarm_residuncer(jj,1) = corr(sort_roibeta, sort_uncer_resid);
    res.corr_roibetarm_absresiduncer(jj,1) = corr(sort_roibeta, abs(sort_uncer_resid));
    res.absresiduncerbytrial{:,jj} = abs(sort_uncer_resid);
    res.logdmabsresiduncer{:,jj} = log(abs(sort_uncer_resid));
    
    btmp = fitglm(zscore(abs(sort_uncer_resid)), zscore(sort_roibeta));
    res.beta_roibetarm_absresiduncer(jj,1) = btmp.Coefficients.Estimate(2);
    btmp = fitglm(zscore(abs(sort_fitresid-mean(sort_fitresid))), zscore(sort_roibeta));
    res.beta_roibetarm_ufituncer(jj,1) = btmp.Coefficients.Estimate(2);

    tmpshuffle = zeros(size(tmpmatrix));
    % shuffle
    nbin = 8;
    bb = round(ntrial/nbin);
    %         binlabel = repmat(1:nbin, bb, 1);
    %         binlabel = binlabel(1:ntrial)';
    for ii = 1:nbin
        if ii~=nbin
            thisidx = bb*(ii-1)+1 : bb*ii;
        else
            thisidx = bb*(ii-1)+1 : ntrial;
        end
        thisbin = tmpmatrix(thisidx, :);
        res.uncerresid_binbyfituncer{jj,ii} = thisbin(:,end-2);
        res.residroibeta_binbyfituncer{jj,ii} = thisbin(:,end-1);
        res.roibeta_residuncer_binbyfituncer{jj,ii} = sortrows([thisbin(:,end-1), abs(thisbin(:,end-2))],2);
        res.ctime_binbyfituncer{jj,ii} = thisbin(:, 5);
        res.fituncer_bybin(jj,ii) = mean(thisbin(:, 6));
        
    end
    
    %% 
    %先把 daccbeta, residualuncer 按fituncer sort、分bin，然后按bin分为edge组
    %和center组，取center组中residualuncer大于70分位数的daccbeta作为centerhigh
    %然后要根据residualuncer的centerlow去找统计上无差异的edgelow
    tmpedge = [res.roibeta_residuncer_binbyfituncer{jj,1}; res.roibeta_residuncer_binbyfituncer{jj,2}; ...
        res.roibeta_residuncer_binbyfituncer{jj,7}; res.roibeta_residuncer_binbyfituncer{jj,8}];
    tmpcenter = [res.roibeta_residuncer_binbyfituncer{jj,3}; res.roibeta_residuncer_binbyfituncer{jj,4}; ...
        res.roibeta_residuncer_binbyfituncer{jj,5}; res.roibeta_residuncer_binbyfituncer{jj,6}];
    
%     subplot(5,6,jj);histogram(tmpedge(:,2));hold on; histogram(tmpcenter(:,2))

%     res.roibeta_edgehigh{jj,1} = tmpedge(tmpedge(:,2)>=prctile(tmpedge(:,2), 60), :);
    initidx = 30;    
    res.roibeta_centrhigh{jj,1} = tmpcenter(tmpcenter(:,2)>=prctile(tmpcenter(:,2), 70), :);
    while initidx
        lowidx = tmpcenter(:,2)>=prctile(tmpcenter(:,2), initidx) & tmpcenter(:,2)<prctile(tmpcenter(:,2), initidx+30); %起始区间 40%-70% 分位数
        res.roibeta_centrlow{jj,1} = tmpcenter(lowidx, :);
        tmpcentrlow = tmpcenter(lowidx, :);
        

        tmpedge(tmpedge(:,2)>max(tmpcentrlow(:,2)),:) = []; %移除tmpedge中所有比centrlow中最大值还大的那些trial
%         tmpedge(tmpedge(:,2)>prctile(tmpcenter(:,2), initidx+20),:) = []; %移除tmpedge中所有比centrlow中最大值的1.5倍高那些trial
%         tmpedge(tmpedge(:,2)>max(tmpcentrlow(:,2))*1.5,:) = [];
        
        nstd = (mean(tmpcentrlow(:,2)) - mean(tmpedge(:,2)))/std(tmpedge(:,2));
        overlaprate = mean(tmpedge(:,2)>=min(tmpcentrlow(:,2)) & tmpedge(:,2)<=max(tmpcentrlow(:,2)));
        if  overlaprate < 0.3 || abs(nstd)>=0.8
            initidx = initidx-1;
        else
            fprintf('%d Overlap rate: %.3f, diststd: %.2f\n',jj, overlaprate, nstd)
            break
        end
    end
    if initidx == 0 || length(tmpedge)< length(tmpcentrlow)
        warning('%d Overlap rate: %.3f\n', jj, overlaprate)
        res.roibeta_edgelow{jj,1} = nan(length(tmpcentrlow),2);
    else
        while true
            tmpidx = randperm(length(tmpedge), length(tmpcentrlow));
            
%             h = ttest2(tmpcentrlow(:,2), tmpedge(tmpidx, 2)); % 因为是不同trial的所以用test2
            h = ttest2(tmpcentrlow(:,2), tmpedge(tmpidx, 2));
            if ~h
                break
            end
        end
        res.roibeta_edgelow{jj,1} = tmpedge(tmpidx, :);
    end
    
    
    %%
    nbin = 4;
    nhf = length(sort_neghf);
    bb = round(nhf/nbin);
    %         binlabel = repmat(1:nbin, bb, 1);
    %         binlabel = binlabel(1:ntrial)';
    for ii = 1:nbin
        if ii~=nbin
            thisidx = bb*(ii-1)+1 : bb*ii;
        else
            thisidx = bb*(ii-1)+1 : nhf;
        end
        thisbin = sort_neghf(thisidx, :);

        res.residroibeta_binbyresiduncer{jj,ii} = thisbin(:,end-1);
        res.residuncer_bybin(jj,ii) = mean(thisbin(:, 7));
    end
    
    nbin = 4;
    nhf = length(sort_poshf);
    bb = round(nhf/nbin);
    %         binlabel = repmat(1:nbin, bb, 1);
    %         binlabel = binlabel(1:ntrial)';
    for ii = 1:nbin
        if ii~=nbin
            thisidx = bb*(ii-1)+1 : bb*ii;
        else
            thisidx = bb*(ii-1)+1 : nhf;
        end
        thisbin = sort_poshf(thisidx, :);

        res.residroibeta_binbyresiduncer{jj,ii+nbin} = thisbin(:,end-1);
        res.residuncer_bybin(jj,ii+nbin) = mean(thisbin(:, 7))*-1;
    end
    
end

getglm = struct2table(getglm);

[res.auccorr1(1,1), res.auccorr1(1,2)] = corr(res.before_aucbeh1, res.before_auc);
[res.auccorr2(1,1), res.auccorr2(1,2)] = corr(res.before_aucbeh2, res.before_auc);

res.uncer_resid_bin_mean = cellfun(@mean, res.uncer_resid_bin);
res.uncer_resid_bin_std = cellfun(@std, res.uncer_resid_bin);
res.uncer_resid_bin_len = cellfun(@length, res.uncer_resid_bin);

res.ctime_binbyrtime_mean = cellfun(@mean, res.ctime_binbyrtime);
res.ctime_binbyrtime_std = cellfun(@std, res.ctime_binbyrtime);
res.ctime_binbyrtime_len = cellfun(@length, res.ctime_binbyrtime);

res.uncerresid_binbyfituncer_mean = cellfun(@mean, res.uncerresid_binbyfituncer);
res.uncerresid_binbyfituncer_std = cellfun(@std, res.uncerresid_binbyfituncer);
res.uncerresid_binbyfituncer_len = cellfun(@length, res.uncerresid_binbyfituncer);
tmpabs = cellfun(@abs, res.uncerresid_binbyfituncer, 'UniformOutput', false);
res.uncerresid_binbyfituncer_absmean = cellfun(@mean, tmpabs);

res.residroibeta_binbyfituncer_mean = cellfun(@mean, res.residroibeta_binbyfituncer);
res.residroibeta_binbyfituncer_std = cellfun(@std, res.residroibeta_binbyfituncer);
res.residroibeta_binbyfituncer_len = cellfun(@length, res.residroibeta_binbyfituncer);

res.ctime_binbyfituncer_mean = cellfun(@mean, res.ctime_binbyfituncer);
res.ctime_binbyfituncer_std = cellfun(@std, res.ctime_binbyfituncer);
res.ctime_binbyfituncer_len = cellfun(@length, res.ctime_binbyfituncer);

res.residroibeta_binbyresiduncer_mean = cellfun(@mean, res.residroibeta_binbyresiduncer);
res.residroibeta_binbyresiduncer_std = cellfun(@std, res.residroibeta_binbyresiduncer);
res.residroibeta_binbyresiduncer_len = cellfun(@length, res.residroibeta_binbyresiduncer);

res.uncer_resid_bin_sem = res.uncer_resid_bin_std ./sqrt(res.uncer_resid_bin_len);

for jj = 1:nsub
    res.corr_bined_meanroibeta_meanuncer(jj,1) = corr(res.residroibeta_binbyresiduncer_mean(jj,:)', ...
                                                    res.residuncer_bybin(jj,:)','type','Spearman');
end

%  res.roibeta_edgehigh_mean = cellfun(@mean,  res.roibeta_edgehigh);
%  res.roibeta_edgelow_mean = cellfun(@mean,  res.roibeta_edgelow);
%  res.roibeta_centrhigh_mean = cellfun(@mean,  res.roibeta_centrhigh);
%  res.roibeta_centrlow_mean = cellfun(@mean,  res.roibeta_centrlow);
    
groupcentrhigh = cell2mat( cellfun(@mean,  res.roibeta_centrhigh,  'UniformOutput',  false));
groupcentrlow = cell2mat( cellfun(@mean,  res.roibeta_centrlow,  'UniformOutput',  false));
groupedgelow = cell2mat( cellfun(@nanmean,  res.roibeta_edgelow,  'UniformOutput',  false));
res.roibeta_mean_split = table( groupcentrhigh(:,1), groupcentrlow(:,1), groupedgelow(:,1), ...
    'VariableNames', {'centrhigh','centrlow','edgelow'});
res.residuncer_mean_split = table( groupcentrhigh(:,2), groupcentrlow(:,2), groupedgelow(:,2), ...
    'VariableNames', {'centrhigh','centrlow','edgelow'});

toc(st)
end


