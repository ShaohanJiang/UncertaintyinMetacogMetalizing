%% AUC correlation across tasks between subjects
clear;home;

load('allBehavTable.mat')

%%
st = tic;
targettask = RDMtbl;
thistask = SAtbl;
maxlen = 10;
ntimes = 1e5;
beforeauc_twotask = [];
for jj = 1:height(RDMtbl)
    thissubno = FindStrinCell(targettask.subnames{jj}, thistask.subnames);
    if length(thissubno)==1
        % target task
        clearidx = targettask.resp{jj} ~=-1 & targettask.conf{jj}~=-1;
        diffc = targettask.diff{jj}(clearidx);
        rtime = targettask.rtime{jj}(clearidx);
        conf = targettask.conf{jj}(clearidx);
        isco = targettask.isc{jj}(clearidx);
        
        dmat = [ones(length(conf),1), diffc, rtime];
        idmat = pinv(dmat);
        beta = idmat*conf;
        conf_resid = conf - dmat*beta;
        
        targetauc = fastAUC(isco*2-1, conf_resid, 1);
        
        % this task
        clearidx = thistask.resp{thissubno} ~=-1 & thistask.conf{thissubno}~=-1;
        diffc = thistask.diff{thissubno}(clearidx);
        rtime = thistask.rtime{thissubno}(clearidx);
        conf = thistask.conf{thissubno}(clearidx);
        isco = thistask.isc{thissubno}(clearidx);
        
        dmat = [ones(length(conf),1), diffc, rtime];
        idmat = pinv(dmat);
        beta = idmat*conf;
        conf_resid = conf - dmat*beta;
        
        thisauc = fastAUC(isco*2-1, conf_resid, 1);
        
        % output
        beforeauc_twotask = [beforeauc_twotask; [thisauc, targetauc]];
        
        % bin & randperm this task
        tmpmatrix = [conf, diffc, rtime, (diffc-mean(diffc)).*(rtime-median(rtime)), isco];
        tmpmatrix = sortrows(tmpmatrix, 4);
        sort_isco = tmpmatrix(:, end);
        ntrial = size(tmpmatrix, 1);
        for bb = 1:maxlen
            tmpshuffle = zeros(size(tmpmatrix));
            % shuffle
            nbin = ceil(ntrial/bb);
            for tt = 1:ntimes
                for ii = 1:nbin
                    if ii~=nbin
                        thisidx = bb*(ii-1)+1 : bb*ii;
                    else
                        thisidx = bb*(ii-1)+1 : ntrial;
                    end
                    thisbin = tmpmatrix(thisidx, :);
                    tmpshuffle(thisidx, :) = thisbin(randperm(length(thisidx)), :);
                end
                shuff_conf = tmpshuffle(:,1);
                shuff_diffrt = tmpshuffle(:,2:3);
                
                dmat = [ones(ntrial,1), shuff_diffrt];
                idmat = pinv(dmat);
                beta = idmat*shuff_conf;
                shuff_conf_resid = shuff_conf - dmat*beta;
                
                afterauc_thistask(jj, bb, tt) = fastAUC(sort_isco*2-1, shuff_conf_resid, 1);
            end
        end
    else
        fprintf('Sub %s didnot perform this task\n', targettask.subnames{jj})
    end
    
end
toc(st)

%%
% erase subs
afterauc_thistask(mean(mean(afterauc_thistask,3),2)==0,:,:) = [];

%
corr_twotask = zeros(ntimes, maxlen);
for tt = 1:ntimes
    corr_twotask(tt,:) = corr(afterauc_thistask(:,:,tt), beforeauc_twotask(:,2))';

    
    
end

%%
figure;
violinplot(corr_twotask, [], 'ViolinColor', [0 0.4470 0.7410], 'ShowData', false, 'ViolinAlpha', 1, ...
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
mbins = mean(corr_twotask);
sembins = std(corr_twotask); %range(CI)/2; %
errorbar(2:10, mbins(2:10), sembins(2:10), '*','Color','k','LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});


%%
meanauc = mean(afterauc_thistask,3);

figure;
violinplot(meanauc, [], 'ViolinColor', [0 0.4470 0.7410], 'ShowData', false, 'ViolinAlpha', 1, ...
    'EdgeColor', [0 0.4470 0.7410], 'BoxColor',  [0 0.4470 0.7410], 'ShowMean', true, 'MedianColor',[0 0.4470 0.7410]);

legend off
box off
set(gca, 'XTick',1:10)
xlim([0.5,10.5])
ylim([0.35, 0.7])
set(gca, 'XTick', 1:10, 'YTick',0.4:0.1:0.7, 'LineWidth',1);
hold on
plot([0.5,10.5], [0.5,0.5], 'k--')

hold on 
mbins = mean(meanauc);
sembins = std(meanauc); %range(CI)/2; %
errorbar(1:10, mbins(1:10), sembins(1:10), '*','Color','k','LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});





