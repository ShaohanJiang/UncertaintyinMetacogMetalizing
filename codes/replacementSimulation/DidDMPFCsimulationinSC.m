%% DMPFC in SC
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
    roibeta_dmpfc_raw = mean(thistbl.DMPFC_JSH_oconf20_raw{jj},2);
    roibeta_dmpfc_rm = mean(thistbl.DMPFC_JSH_oconf20_rmDiffRT{jj},2);
    
    if length(uncer) == length(roibeta_dmpfc_rm)
        ntrial = length(uncer);
    end
    % Before
    [~,~,stats] = glmfit([diffc, rtime], uncer);
    before_corr(jj,1) = corr(stats.resid, roibeta_dmpfc_rm);
    before_auc(jj,1) = fastAUC(isco*2-1, roibeta_dmpfc_rm, -1);
    
    tmpmatrix = [uncer, diffc, rtime, (diffc-mean(diffc)).*(rtime-median(rtime)), roibeta_dmpfc_rm, isco];
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
%%

save res_dmpfcinsc_simulation_replacement_10k after* before* maxlen

%%

after_corrauc = zeros(ntimes, maxlen);

    for jj = 1:ntimes
        after_corrauc(jj, :) = corr(before_auc, after_auc(:,:,jj));

    end

figure; histogram(after_corrauc(:,10))
box off
hold on; plot([0,0], [0,700], '--', 'LineWidth',1.5)

%% Painting

bncorr_bins = mean(after_corr, 3);


mbins = mean(bncorr_bins);
sembins = std(bncorr_bins)/sqrt(30);
figure;
errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black')
% set(gca,'XTick',1:4, 'YTick',0.5:0.1:0.9, 'LineWidth',2.5, 'XTickLabel',{}, 'YTickLabel',{})
box off
set(gca, 'XTick', 1:10, 'YTick',0:0.02:0.1, 'LineWidth',1)
xlim([0.5,10.5])
set(gca, 'XTickLabel',{}, 'YTickLabel',{});

%%
neauc_bins = mean(after_auc, 3);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(30);
figure;
errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black')
% set(gca,'XTick',1:4, 'YTick',0.5:0.1:0.9, 'LineWidth',2.5, 'XTickLabel',{}, 'YTickLabel',{})
box off
xlim([0.5,10.5])
ylim([0.49, 0.54])
set(gca, 'XTick', 1:10, 'YTick',0.48:0.01:0.54, 'LineWidth',1);
hold on
plot([0.5, 10.5], [0.5, 0.5], 'k--')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});


%%
after_auc(:,5,:) = tmp.after_auc(:,5,:);
after_corr(:,5,:) = tmp.after_corr(:,5,:);
after_corr(:,4,:) = tmp.after_corr(:,4,:);
after_corr(:,3,:) = tmp.after_corr(:,3,:);
after_corr(:,2,:) = tmp.after_corr(:,2,:);
after_corr(:,1,:) = tmp.after_corr(:,1,:);


%%
figure;
violin(after_corrauc, 'facecolor',[0 0.4470 0.7410],'edgecolor','none',...  % 1 0.56 0
'bw',0.3,'mc','k','medc',[], 'facealpha',0.75);
legend off
box off
set(gca, 'XTick',1:10)
xlim([0.5,10.5])
ylim([-1.75, 0.54])
set(gca, 'XTick', 1:10, 'YTick',0.48:0.01:0.54, 'LineWidth',2.5);
hold on

plot([0,11], [0,0], 'k--')

%%
for ii = 1:size(after_corrauc,2)
    tmp = randn(size(after_corrauc, 1),1)*std(after_corrauc(:,ii));
    p(1,ii) = kstest2(after_corrauc(:,ii), tmp);
    [~, p(2,ii), CI(ii)] = ztest(after_corrauc(:,ii), 0, std(after_corrauc(:,ii)));
    
end
%%
figure;
violinplot(after_corrauc, [], 'ViolinColor', [0 0.4470 0.7410], 'ShowData', false, 'ViolinAlpha', 1, ...
    'EdgeColor', [0 0.4470 0.7410], 'BoxColor',  [0 0.4470 0.7410], 'ShowMean', true, 'MedianColor',[0 0.4470 0.7410]);

legend off
box off
set(gca, 'XTick',1:10)
xlim([0.5,10.5])
ylim([-1, 1])
set(gca, 'XTick', 1:10, 'YTick',-1:0.5:1, 'LineWidth',1);
hold on
plot([0.5,10.5], [0,0], 'k--')

hold on 
mbins = mean(after_corrauc);
sembins =  range(CI)/2; %std(after_corrauc)
errorbar(2:10, mbins(2:10), sembins(2:10), '.','Color','k','LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});
