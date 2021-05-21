%% PAINT RESULT

%% ROI Decoding correctness
cd E:\Server_Backup\deconv20
load('res_decoding_correctness.mat')

accu_rdm = mean(correctness_rdm_4roi_pca2, 3)-0.5;
accu_oc = mean(correctness_oc_4roi_pca2, 3)-0.5;
% accu_rdm(:,3) = [];
% accu_oc(:,3) = [];

meanall=[mean(accu_rdm); mean(accu_oc)];
semall=[std(accu_rdm)./sqrt(28); std(accu_oc)./sqrt(30)];

colors = [0,0,1; 1 0 0];
[bar_xtick, hb, he] = errorbar_groups(meanall, semall, 'bar_colors', colors, 'bar_width', 0.4, 'bar_interval', 0.08,...
    'optional_bar_arguments',{'LineWidth',1.5},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',4.5, 'LineWidth',1.4});

hold on
% RDM
b = get(hb(1));
x = repmat(b.XData, length(accu_rdm), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), accu_rdm(:,jj), 'k.')
end
% OC
b = get(hb(2));
x = repmat(b.XData, length(accu_oc), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), accu_oc(:,jj), 'k.')
end
b.BaseLine.LineWidth = 1.5;

hold off

ylim([-0.23, 0.2])
set(gca, 'YTick',-0.2:0.1:0.2, 'YTickLabel',{},'XTickLabel',{}, 'LineWidth',1.5)

%% neutral ROI Decoding correctness
cd E:\Server_Backup\deconv20
load('res_correctness_rdmoc_neurosynth.mat')

clear accu_rdm accu_oc
accu_rdm_allroi = mean(correctness_rdm_4roi_pcatmp51, 3)-0.5;
accu_oc_allroi = mean(correctness_oc_4roi_pcatmp51, 3)-0.5;
accu_rdm = [accu_rdm_allroi(:,1), accu_rdm_allroi(:,3), accu_rdm_allroi(:,4), accu_rdm_allroi(:,5)];
accu_oc = [accu_oc_allroi(:,1), accu_oc_allroi(:,3), accu_oc_allroi(:,4), accu_oc_allroi(:,5)];

meanall=[mean(accu_rdm); mean(accu_oc)];
semall=[std(accu_rdm)./sqrt(size(accu_rdm,1)); std(accu_oc)./sqrt(size(accu_oc,1))];

colors = [0,0,1; 1 0 0];
[bar_xtick, hb, he] = errorbar_groups(meanall, semall, 'bar_colors', colors, 'bar_width', 0.4, 'bar_interval', 0.08,...
    'optional_bar_arguments',{'LineWidth',1.5},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',4.5, 'LineWidth',1.4});

hold on
% RDM
b = get(hb(1));
x = repmat(b.XData, length(accu_rdm), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), accu_rdm(:,jj), 'k.')
end
% OC
b = get(hb(2));
x = repmat(b.XData, length(accu_oc), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), accu_oc(:,jj), 'k.')
end
b.BaseLine.LineWidth = 1.5;

hold off

ylim([-0.23, 0.2])
set(gca, 'YTick',-0.2:0.1:0.2, 'YTickLabel',{},'XTickLabel',{}, 'LineWidth',1.5)

%% ROI AUC
cd E:\Server_Backup\deconv20
load('res_roiauc_ourroi_resid.mat')

meanall=[mean(roiaucresid_rdm); mean(roiaucresid_oc)]-0.5;
semall=[std(roiaucresid_rdm)./sqrt(28); std(roiaucresid_oc)./sqrt(30)];

colors = [0,0,1; 1 0 0];
[bar_xtick, hb, he] = errorbar_groups(meanall, semall, 'bar_colors', colors, 'bar_width', 0.4, 'bar_interval', 0.08,...
    'optional_bar_arguments',{'LineWidth',1.5},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',4.5, 'LineWidth',1.4});

hold on
% RDM
b = get(hb(1));
x = repmat(b.XData, length(roiaucresid_rdm), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), roiaucresid_rdm(:,jj)-0.5, 'k.')
end
% OC
b = get(hb(2));
x = repmat(b.XData, length(roiaucresid_oc), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), roiaucresid_oc(:,jj)-0.5, 'k.')
end
b.BaseLine.LineWidth = 1.5;
hold off

ylim([-0.23,0.2])
set(gca, 'YTick',-0.2:0.1:0.2, 'YTickLabel',{},'XTickLabel',{}, 'LineWidth',1.5)

%% ROI AUC (Neutral ROI)
cd E:\Server_Backup\deconv20
load('res_roiauc_neurosynth_rdmoc.mat')

roiaucresid_rdm = roiauc_rdm_resid(:,[1,2,3,4]);
roiaucresid_oc = roiauc_oc_resid(:,[1,2,3,4]);
meanall=[mean(roiaucresid_rdm); mean(roiaucresid_oc)]-0.5;
semall=[std(roiaucresid_rdm)./sqrt(28); std(roiaucresid_oc)./sqrt(30)];

colors = [0,0,1; 1 0 0];
[bar_xtick, hb, he] = errorbar_groups(meanall, semall, 'bar_colors', colors, 'bar_width', 0.4, 'bar_interval', 0.08,...
    'optional_bar_arguments',{'LineWidth',1.5},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',4.5, 'LineWidth',1.4});

hold on
% RDM
b = get(hb(1));
x = repmat(b.XData, length(roiaucresid_rdm), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), roiaucresid_rdm(:,jj)-0.5, 'k.')
end
% OC
b = get(hb(2));
x = repmat(b.XData, length(roiaucresid_oc), 1);
for jj = 1:size(x,2)
    plot(x(:,jj), roiaucresid_oc(:,jj)-0.5, 'k.')
end
b.BaseLine.LineWidth = 1.5;
hold off

ylim([-0.23,0.2])
set(gca, 'YTick',-0.2:0.1:0.2, 'YTickLabel',{},'XTickLabel',{}, 'LineWidth',1.5)

%% Plot CTIME across conf
cd E:\Server_Backup\MM_behavedata
load('allBehavTable.mat')

sem = @(x) std(x)/sqrt(size(x,1));

% [OC; SC; OA; SA]
%OC
thistbl = OCtbl;
allmedct = [];
for jj = 1:height(thistbl)
    clearidx = thistbl.resp{jj} ~=-1 & thistbl.conf{jj}~=-1;
    thisconf = thistbl.conf{jj}(clearidx);
    thisct = thistbl.ctime{jj}(clearidx);
    allmedct(jj,:) = [median(thisct(thisconf==1)), median(thisct(thisconf==2)), median(thisct(thisconf==3)), median(thisct(thisconf==4))];

end
meanforplot(1,:) = mean(allmedct);
semforplot(1,:) = sem(allmedct);
% SC
thistbl = SCtbl;
allmedct = [];
for jj = 1:height(thistbl)
    clearidx = thistbl.resp{jj} ~=-1 & thistbl.conf{jj}~=-1;
    thisconf = thistbl.conf{jj}(clearidx);
    thisct = thistbl.ctime{jj}(clearidx);
    allmedct(jj,:) = [median(thisct(thisconf==1)), median(thisct(thisconf==2)), median(thisct(thisconf==3)), median(thisct(thisconf==4))];

end
meanforplot(2,:) = mean(allmedct);
semforplot(2,:) = sem(allmedct);
% OA
thistbl = OAtbl;
allmedct = [];
for jj = 1:height(thistbl)
    clearidx = thistbl.resp{jj} ~=-1 & thistbl.conf{jj}~=-1;
    thisconf = thistbl.conf{jj}(clearidx);
    thisct = thistbl.ctime{jj}(clearidx);
    allmedct(jj,:) = [median(thisct(thisconf==1)), median(thisct(thisconf==2)), median(thisct(thisconf==3)), median(thisct(thisconf==4))];

end
meanforplot(3,:) = mean(allmedct);
semforplot(3,:) = sem(allmedct);
% SA
thistbl = SAtbl;
allmedct = [];
for jj = 1:height(thistbl)
    clearidx = thistbl.resp{jj} ~=-1 & thistbl.conf{jj}~=-1;
    thisconf = thistbl.conf{jj}(clearidx);
    thisct = thistbl.ctime{jj}(clearidx);
    allmedct(jj,:) = [median(thisct(thisconf==1)), median(thisct(thisconf==2)), median(thisct(thisconf==3)), median(thisct(thisconf==4))];

end
meanforplot(4,:) = mean(allmedct);
semforplot(4,:) = sem(allmedct);


figure;
taskcolors = [
    1 0 0; ... OC
    1 0 1; ... SC
    0 1 0; ... OA
    0,1,1]; % SA
for ii = 1:size(meanforplot,1)

    eb = errorbar([1,2,3,4], meanforplot(ii,:), semforplot(ii,:), 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
    eb.Color = taskcolors(ii,:);
    hold on
    % title('OA')
end
box off
set(gca,'XTick',1:4, 'YTick',0.5:0.1:0.9, 'LineWidth',2.5, 'XTickLabel',{}, 'YTickLabel',{})
xlim([0.8,4.2])
ylim([0.5,0.9])

%%
cd E:\Server_Backup\MM_meants\
load('taskroibeta.mat')

thisdata = taskroibeta_IPLrtime6mm; % taskroibeta_IPLrtime6mm; taskroibeta_IPLdiff6mm

allmean = cellfun(@mean, thisdata);
allsem = cellfun(@std, thisdata)./sqrt(cellfun(@length, thisdata));

allmean(:,1) = [];
allsem(:,1) = [];
taskcolors = [
    1 0 0; ... OC
    1 0 1; ... SC
    0 1 0; ... OA
    0,1,1]; % SA

figure;
b = bar(allmean, 'LineWidth',2, 'BarWidth',0.35);
b.FaceColor = 'flat';
b.CData = taskcolors;
hold on
eb = errorbar([1,2,3,4], allmean, allsem, 'LineStyle', 'none','LineWidth', 2, 'Marker','o', 'MarkerFaceColor','white','MarkerEdgeColor','black');
eb.Color = 'black';
box off
ax = get(gca);
set(gca,'XTick',1:4, 'YTick', min(ax.YTick):max(ax.YTick)/4:max(ax.YTick), 'LineWidth',2.5, 'XTickLabel',{}, 'YTickLabel',{})


%% ttest correction
cd E:\Server_Backup\deconv20
% roi from 
load('res_roiauc_ourroi_resid.mat')
load('res_decoding_correctness.mat')
p_bias = table();
% decoding
accu_rdm = mean(correctness_rdm_4roi_pca2, 3)-0.5;
accu_oc = mean(correctness_oc_4roi_pca2, 3)-0.5;
accu_rdm(:,3) = [];
accu_oc(:,3) = [];
[~,p] = ttest(accu_rdm,0,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_bias.twotailrdmdecoding = p';
p_bias.twotailrdmdecoding_fdr = fdr';
[~,p] = ttest(accu_oc,0,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_bias.twotailocdecoding = p';
p_bias.twotailocdecoding_fdr = fdr';
[~,p] = ttest(accu_rdm,0,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_bias.onetailrdmdecoding = p';
p_bias.onetailrdmdecoding_fdr = fdr';
[~,p] = ttest(accu_oc,0,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_bias.onetailocdecoding = p';
p_bias.onetailocdecoding_fdr = fdr';
% roiauc
[~,p] = ttest(roiaucresid_rdm, 0.5,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_bias.twotailrdmroiauc = p';
p_bias.twotailrdmroiauc_fdr = fdr';
[~,p] = ttest(roiaucresid_oc,0.5,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_bias.twotailocroiauc = p';
p_bias.twotailocroiauc_fdr = fdr';
[~,p] = ttest(roiaucresid_rdm, 0.5,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_bias.onetailrdmroiauc = p';
p_bias.onetailrdmroiauc_fdr = fdr';
[~,p] = ttest(roiaucresid_oc,0.5,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_bias.onetailocroiauc = p';
p_bias.onetailocroiauc_fdr = fdr';

% roi from neurosynth
load('res_correctness_rdmoc_neurosynth.mat')
load('res_roiauc_neurosynth_rdmoc.mat')
clear accu_rdm accu_oc roiaucresid_rdm roiaucresid_oc
p_neurosynth = table();
% decoding
accu_rdm_allroi = mean(correctness_rdm_4roi_pcatmp51, 3)-0.5;
accu_oc_allroi = mean(correctness_oc_4roi_pcatmp51, 3)-0.5;
accu_rdm = [accu_rdm_allroi(:,1), accu_rdm_allroi(:,3), accu_rdm_allroi(:,4), accu_rdm_allroi(:,5)];
accu_oc = [accu_oc_allroi(:,1), accu_oc_allroi(:,3), accu_oc_allroi(:,4), accu_oc_allroi(:,5)];
[~,p] = ttest(accu_rdm,0,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.twotailrdmdecoding = p';
p_neurosynth.twotailrdmdecoding_fdr = fdr';
[~,p] = ttest(accu_oc,0,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.twotailocdecoding = p';
p_neurosynth.twotailocdecoding_fdr = fdr';
[~,p] = ttest(accu_rdm,0,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.onetailrdmdecoding = p';
p_neurosynth.onetailrdmdecoding_fdr = fdr';
[~,p] = ttest(accu_oc,0,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.onetailocdecoding = p';
p_neurosynth.onetailocdecoding_fdr = fdr';
% roiauc
roiaucresid_rdm = roiauc_rdm_resid(:,[1,2,3,4]);
roiaucresid_oc = roiauc_oc_resid(:,[1,2,3,4]);
[~,p] = ttest(roiaucresid_rdm, 0.5,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.twotailrdmroiauc = p';
p_neurosynth.twotailrdmroiauc_fdr = fdr';
[~,p] = ttest(roiaucresid_oc,0.5,'Tail','both');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.twotailocroiauc = p';
p_neurosynth.twotailocroiauc_fdr = fdr';
[~,p] = ttest(roiaucresid_rdm, 0.5,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.onetailrdmroiauc = p';
p_neurosynth.onetailrdmroiauc_fdr = fdr';
[~,p] = ttest(roiaucresid_oc,0.5,'Tail','right');
fdr = mafdr(p, 'BHFDR',true);
p_neurosynth.onetailocroiauc = p';
p_neurosynth.onetailocroiauc_fdr = fdr';
% output
 writetable(p_bias, 'p_bias.csv','FileType','text')
  writetable(p_neurosynth, 'p_neurosynth.csv','FileType','text')
  
%%
  