%%
load('figS2_data.mat')
all_res = {resrdm, resoc, ressc, resoa, ressa};
taskcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];

%% plot mrtime maccu split by difficulty

for ii = 1:1
    thisdata = all_res{ii}.maccu*100;
    matmean = mean(thisdata);
    matsem = std(thisdata)./sqrt(length(thisdata));
    fig = figure;
    set(gcf,'unit','centimeters','position',[3 5 5 4])
    yyaxis left
    e = errorbar(matmean, matsem, '-','LineWidth',1,'MarkerSize',0.1,'MarkerFaceColor','white','CapSize',2);
    set(gca, 'linewidth', 1, 'YTick',25:25:100);
    ylim([25, 100])
    hold on
    
    thisdata = all_res{ii}.mrtall;
    matmean = mean(thisdata);
    matsem = std(thisdata)./sqrt(length(thisdata));
%     fig = figure;
%     set(gcf,'unit','centimeters','position',[3 5 4.52 3.45])
    yyaxis right
    e = errorbar(matmean, matsem, '-','LineWidth',1,'MarkerSize',0.1,'MarkerFaceColor','white','CapSize',2);
    set(gca, 'linewidth', 1, 'YTick',0.5:0.1:0.9);
    ylim([0.5, 0.9])
    xlim([0.5, 4.5])
    box off
end


%% plot mconf split by difficulty and isc
rtcolor = [166,54,3; 253,190,133]./255;
for ii = 1:1
    thisdata = all_res{ii}.mconf;
    ctmp = 5-squeeze(thisdata(:,1,:))'; % uncertainty
    etmp = 5-squeeze(thisdata(:,2,:))'; % uncertainty
    meanc = nanmean(ctmp);
    meane = nanmean(etmp);
    semc = nanstd(ctmp)./sqrt(sum(~isnan(ctmp)));
    seme = nanstd(etmp)./sqrt(sum(~isnan(etmp)));
    
    matmean = [meanc; meane];
    matsem = [semc; seme];
    fig = figure;
    set(gcf,'unit','centimeters','position',[3 5 5 4])
    for jj = 1:2
        e = errorbar(matmean(jj,:), matsem(jj,:), '-','LineWidth',1,'MarkerSize',0.1,'MarkerFaceColor','white','CapSize',2);
        % set(gca, 'linewidth', 1, 'YTick',-0.2:0.2:0.5);
        e.Color = rtcolor(jj,:);
        hold on
    end
    ylim([1, 4])
    xlim([0.5, 4.5])
    box off
end

%% plot mrtime split by difficulty and isc

rtcolor = [166,54,3; 253,190,133]./255;
for ii = 1:1
    thisdata = all_res{ii}.mrtime;
    ctmp = squeeze(thisdata(:,1,:))';
    etmp = squeeze(thisdata(:,2,:))';
    meanc = nanmean(ctmp);
    meane = nanmean(etmp);
    semc = nanstd(ctmp)./sqrt(sum(~isnan(ctmp)));
    seme = nanstd(etmp)./sqrt(sum(~isnan(etmp)));
    
    matmean = [meanc; meane];
    matsem = [semc; seme];
    fig = figure;
    set(gcf,'unit','centimeters','position',[3 5 5 4])
    for jj = 1:2
        e = errorbar(matmean(jj,:), matsem(jj,:), '-','LineWidth',1,'MarkerSize',0.1,'MarkerFaceColor','white','CapSize',2);
        e.Color = rtcolor(jj,:);
        hold on
    end
    set(gca, 'YTick', 0.4:0.2:1.2);
    ylim([0.4, 1.2])
    xlim([0.5, 4.5])
    box off
end
