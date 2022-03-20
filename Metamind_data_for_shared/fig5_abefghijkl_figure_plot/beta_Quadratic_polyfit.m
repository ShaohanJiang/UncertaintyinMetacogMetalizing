%%
load('fig5metadata.mat')

all_res = {resrdm, resoc, ressc, resoa, ressa};
taskcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];


%% fit for fig 5a, binned ctime test
clear allP
for ii = 1:5
    tmpa = all_res{ii}.ctime_bin;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end
matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-1.2:0.4:0);
box off
ylim([-1.2, 0])


% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})

%% fit for fig 5b ctime_binbyfituncer_mean


for ii = 1:5
    tmpa = all_res{ii}.ctime_binbyfituncer_mean;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end

matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-1.2:0.4:0);
box off
ylim([-1.2, 0])

% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})


%% fit for fig 5d PLOT binned absuncer (revise) 


for ii = 1:5
    tmpa = all_res{ii}.absuncer_binbyfituncer_mean;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end

% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})

matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-1.2:0.4:0);
box off
ylim([-1.2, 0])

%% fit for fig 5g

for ii = 1:5
    tmpa = all_res{ii}.uncerresid_binbyfituncer_std;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end

% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})

matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-1.2:0.4:0);
box off
ylim([-1.2, 0])

%% fit for fig 5h
for ii = 1:5
    tmpa = all_res{ii}.residroibeta_binbyfituncer_mean;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end

% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})

matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-1.2:0.4:0);
box off
ylim([-1.2, 0])

%% fit for fig 5j
for ii = 1:5
    tmpa = all_res{ii}.uncerresid_binbyfituncer_absmean;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end

% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})

matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',-1.2:0.4:0);
box off
ylim([-1.2, 0])

%% fit for fig 5k

for ii = 1:5
    tmpa = all_res{ii}.residroibeta_binbyresiduncer_mean;
    P = Get_QuaPolyfit(tmpa);
    allP{ii} = P(:,1);
end

% [~,p,~,stat] = ttest(allP{1})
% [~,p,~,stat] = ttest(allP{2})
% [~,p,~,stat] = ttest(allP{3})
% [~,p,~,stat] = ttest(allP{4})
% [~,p,~,stat] = ttest(allP{5})

matmean = cellfun(@mean, allP);
matsem = bsxfun(@rdivide, cellfun(@std, allP), sqrt(cellfun(@length, allP)));

fig = figure;
set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
[bar_xtick, hb, he] = errorbar_groups(matmean, matsem, 'bar_width', 0.7, 'bar_interval', 0.08, 'FigID', fig, ...
    'optional_bar_arguments',{'LineWidth',1, 'FaceColor','flat', 'CData',taskcolor},...
    'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'MarkerSize', 3,'CapSize',3, 'LineWidth',1});
b = get(hb);
b.BaseLine.LineWidth = 1.2;
set(gca, 'linewidth', 1, 'YTick',0:0.2:0.6, 'XTickLabel', {});
box off
ylim([0, 0.6])





