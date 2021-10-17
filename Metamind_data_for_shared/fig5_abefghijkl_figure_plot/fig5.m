%%
load('fig5metadata.mat')
all_res = {resrdm, resoc, ressc, resoa, ressa};
taskcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];

%% fig5a ctime_bin

neauc_bins = all_res{2}.ctime_bin;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'red';

neauc_bins = all_res{3}.ctime_bin;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'magenta';

neauc_bins = all_res{4}.ctime_bin;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = 'green';

neauc_bins = all_res{5}.ctime_bin;
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
hold on;
e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
e.Color = [0 1 1];

set(gca,'XTick',1:4, 'YTick',0.5:0.1:1, 'LineWidth',2.5)
box off
xlim([0.5,4.5])
ylim([0.5,1])
title('ctime binned by uncertainty')

%% fig 5b ctime_binbyfituncer_mean

figure; 
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.ctime_binbyfituncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

set(gca,'XTick',1:nbin, 'YTick',0.5:0.1:1, 'LineWidth',1.5)
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
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.residroibeta_binbyfituncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

% set(gca,'XTick',1:4, 'LineWidth',1.5)
set(gca,'XTick',1:nbin, 'YTick',-60:20:60,  'LineWidth',1.4)
box off
xlim([0.5,nbin+0.5])
ylim([-60, 40])

%% ���� fig 5g beta_roibetarm_ufituncer

for jj = 1:length(all_res)
    thisr{jj} = all_res{jj}.beta_roibetarm_ufituncer;
end

matmean = cellfun(@mean, thisr);
matsem = bsxfun(@rdivide, cellfun(@std, thisr), sqrt(cellfun(@length, thisr)));

[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.4, 'bar_interval', 0.08,...
    'optional_bar_arguments',{'LineWidth',1.5, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',4.5, 'LineWidth',1.4});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1.2, 'YTick',-0.2:0.1:0.2);
box off
ylim([-0.2, 0])

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
for ii = 1:length(all_res)
    thisres = all_res{ii};
    neauc_bins = thisres.residroibeta_binbyresiduncer_mean;
    [nsub, nbin] = size(neauc_bins);
    mbins = mean(neauc_bins);
    sembins = std(neauc_bins)/sqrt(nsub);
    e = errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black');
    e.Color = taskcolor(ii,:);
    hold on
end
hold off

% set(gca,'XTick',1:4, 'LineWidth',1.5)
set(gca,'XTick',1:nbin, 'YTick',-60:20:60,  'LineWidth',1.4)
box off
xlim([0.5,nbin+0.5])
ylim([-60, 40])

%% ���� fig 5j beta_roibetarm_absresiduncer

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
