%%
load('roibetatsindifferentpara.mat')
 
nroi = width(Tdiff);
sem = @(x) std(x)./sqrt(length(x));
colors = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];
clear atmp btmp ctmp
for ii = 1:nroi
    tmproibeta = table2array([Tdiff(:,ii), Trtime(:,ii), Tuncer(:,ii)]);
    
    
    for jj = 1:3
%         roistats = table;
        atmp = splitroibeta(tmproibeta(:,jj), tasknsub2);
        ctmp = atmp;
        atmp = {ctmp{1}, ctmp{3}, ctmp{2}, ctmp{5}, ctmp{4}};
        for mm = 1:5
            btmp{mm} = atmp{mm}/(std(atmp{mm})./sqrt(length(atmp{mm})));
%             [~,p,~,stats] = ttest(btmp{mm});
%             roistats.p(mm) = p;
%             roistats.t(mm) = stats.tstat;
%             roistats.df(mm) = stats.df;
        end
        
        meanmat(jj,:) = cellfun(@mean, btmp);
        semmat(jj,:) = cellfun(@(x) std(x)./sqrt(length(x)), btmp);
    end
    
    
    [bar_xtick, hb, he] = errorbar_groups(meanmat', semmat', 'bar_colors', colors, 'bar_width', 0.4, 'bar_interval', 0.25,...
        'optional_bar_arguments',{'LineWidth',1.5},...
        'optional_errorbar_arguments',{'LineStyle','none', 'Marker','o','MarkerFaceColor','white', 'CapSize',1, 'LineWidth',1.4});
    
%     names = split(Tdiff.Properties.VariableNames{ii}, '_');
%     title(names{3})
    ylim([-4, 8])
    set(gca, 'ytick', -4:2:8, 'XTickLabel', {}, 'yticklabel', {}, 'linewidth', 1.3)
    b = get(hb(1));
    b.BaseLine.LineWidth = 1.3;
end