%%
load('fig5metadata.mat')
all_res = {resrdm, resoc, ressc, resoa, ressa};
taskcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];

%% fig5a ctime_bin
figure;
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
for ii = 2:5
neauc_bins = all_res{ii}.ctime_bin;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',1,'MarkerSize',3,'MarkerFaceColor','white','MarkerEdgeColor',taskcolor(ii,:),'CapSize',2);
e.Color = taskcolor(ii,:);
end

set(gca,'XTick',1:4, 'YTick',0.5:0.1:1, 'LineWidth',1)
box off
xlim([0.75,4.25])
ylim([0.5,1])
% title('ctime binned by uncertainty')
set(gca, 'XTickLabel', {}, 'YTickLabel',{});

%% fig 5b ctime_binbyfituncer_mean

figure;
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.ctime_binbyfituncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',1,'MarkerSize',3,'MarkerFaceColor','white','MarkerEdgeColor',taskcolor(ii,:),'CapSize',2);
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

set(gca,'XTick',1:nbin, 'YTick',0.5:0.1:1, 'LineWidth',1)
box off
xlim([0.5,nbin+0.5])
ylim([0.5, 1])

%% fig 5e uncerresid_binbyfituncer_std
neauc_bins = all_res{1}.uncerresid_binbyfituncer_std;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
figure;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'blue';

neauc_bins = all_res{2}.uncerresid_binbyfituncer_std;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'red';

neauc_bins = all_res{3}.uncerresid_binbyfituncer_std;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'magenta';

neauc_bins = all_res{4}.uncerresid_binbyfituncer_std;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'green';

neauc_bins = all_res{5}.uncerresid_binbyfituncer_std;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = [0 1 1];

set(gca,'XTick',1:8, 'YTick',0:0.25:1, 'LineWidth',1.5)
box off
xlim([0.5, 8.5])
ylim([0, 1])

%% fig 5f residroibeta_binbyfituncer_mean

figure; 
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.residroibeta_binbyfituncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',1,'MarkerSize',3,'MarkerFaceColor','white','MarkerEdgeColor',taskcolor(ii,:),'CapSize',2);
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

% set(gca,'XTick',1:4, 'LineWidth',1.5)
set(gca,'XTick',1:nbin, 'YTick',-70:35:70,  'LineWidth',1)
box off
xlim([0.5,nbin+0.5])
ylim([-70, 35])
set(gca, 'XTickLabel', {}, 'YTickLabel',{});

%% fig 5g beta_resid_meanroibeta_stduncer_bined

for jj = 1:length(all_res)
    thisr{jj} = all_res{jj}.beta_resid_meanroibeta_stduncer_bined;
%     [~, p, ~, t] = ttest(thisr{jj})
end

matmean = cellfun(@mean, thisr);
matsem = bsxfun(@rdivide, cellfun(@std, thisr), sqrt(cellfun(@length, thisr)));

fig = figure;
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-0.2:0.2:0.5);
box off
ylim([0, 0.5])

% fig = figure;
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
% xticks = 0.85+(0:4).*1.2;
% b= bar(xticks, matmean, 'BarWidth', 0.62);
% b.FaceColor = 'flat';
% b.CData = taskcolor;
% 
% hold on
% % xlabel(rowname')
% errorbar(xticks,matmean, matsem, 'linestyle', 'none','LineWidth',1,'Color','k', ...
%     'CapSize',2, 'MarkerEdgeColor','black','Marker','o', 'MarkerSize',3, 'MarkerFaceColor','white'); %
% hold off
% 
% set(gca,'ytick', 0:0.1:0.5);
% ylim([0, 0.5])
% xlim([0, 6.6])
% box off
% set(gca, 'XTickLabel', {}, 'YTickLabel',{});

%% fig 5h uncerresid_binbyfituncer_absmean

figure;
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.uncerresid_binbyfituncer_absmean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

set(gca,'XTick',1:nbin, 'YTick',0:0.2:1, 'LineWidth',1.5)
box off
xlim([0.5,nbin+0.5])
ylim([0, 1])

%% fig 5i residroibeta_binbyresiduncer_mean

figure; 
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.residroibeta_binbyresiduncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',1,'MarkerSize',3,'MarkerFaceColor','white','MarkerEdgeColor',taskcolor(ii,:),'CapSize',2);
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

% set(gca,'XTick',1:4, 'LineWidth',1.5)
set(gca,'XTick',1:nbin, 'YTick',-50:25:50,  'LineWidth',1)
box off
xlim([0.5,nbin+0.5])
ylim([-50, 25])

%% бшнд fig 5j beta_roibetarm_absresiduncer

for jj = 1:length(all_res)
    thisr{jj} = all_res{jj}.beta_roibetarm_absresiduncer;
end

matmean = cellfun(@mean, thisr);
matsem = bsxfun(@rdivide, cellfun(@std, thisr), sqrt(cellfun(@length, thisr)));

[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.4, 'bar_interval', 0.08,...
    'optional_bar_arguments',{'LineWidth',1.5, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',4.5, 'LineWidth',1.4});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
% xlabel('RDM', 'OC', 'SC', 'OA', 'SA')
% set(gca,'XTickLabel',{'RDM', 'OC', 'SC', 'OA', 'SA'})
set(gca, 'linewidth', 1.2, 'YTick',-0.2:0.1:0.2);
box off
ylim([0,0.2])


%% fig 5k & L roibeta_mean_split, residuncer_mean_split
tasktitle = {'RDM', 'OC', 'SC', 'OA', 'SA'};
taskcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];

for jj = 1:length(all_res)
    tmp = table2array(all_res{jj}.roibeta_mean_split);
    tmp(isnan(tmp(:,3)),:) = [];
    nsub = length(tmp);
    figure; PlotErrorbar(tmp)
    title(tasktitle{jj})
    xticklabels(all_res{jj}.roibeta_mean_split.Properties.VariableNames)
    fprintf('-------------\n')
%     ylim([-50,30])
end

for jj = 1:length(all_res)
    tmp = table2array(all_res{jj}.residuncer_mean_split);
    tmp(isnan(tmp(:,3)),:) = [];
    nsub = length(tmp);
    figure; PlotErrorbar(tmp)
    title(tasktitle{jj})
    xticklabels(all_res{jj}.roibeta_mean_split.Properties.VariableNames)
    fprintf('-------------\n')
%     ylim([-50,30])
end

%% PLOT binned absuncer (revise)
% all_res = {res_rdm_dacc, res_oc_dacc, res_sc_dacc, res_oa_dacc, res_sa_dacc};
taskcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];

figure;
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
for ii = 2:5
neauc_bins = all_res{ii}.absuncer_bin;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',1,'MarkerSize',3,'MarkerFaceColor','white','MarkerEdgeColor',taskcolor(ii,:),'CapSize',2);
e.Color = taskcolor(ii,:);
end

set(gca,'XTick',1:4, 'YTick',-2:1:0, 'LineWidth',1)
box off
xlim([0.75,4.25])
ylim([-2,0])
% title('absuncer binned by uncertainty')

figure;
% set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.absuncer_binbyfituncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',1,'MarkerSize',3,'MarkerFaceColor','white','MarkerEdgeColor',taskcolor(ii,:),'CapSize',2);
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

set(gca,'XTick',1:nbin, 'YTick',-2:1:0, 'LineWidth',1)
box off
xlim([0.5,nbin+0.5])
ylim([-2, -0])
