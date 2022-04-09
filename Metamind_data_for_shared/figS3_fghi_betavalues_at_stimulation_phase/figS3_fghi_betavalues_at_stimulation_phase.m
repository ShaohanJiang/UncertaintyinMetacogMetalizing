%%
load('roibetaindifferpara_osti.mat')
 
tasknsub2 = [28, 30, 30, 28, 28];
nroi = width(Tostiuncer);
sem = @(x) std(x)./sqrt(length(x));
colors = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];
clear atmp btmp ctmp
roistats = table;
for ii = 1:nroi
    tmproibeta = table2array(Tostiuncer(:,ii));
    names = split(Tostiuncer.Properties.VariableNames{ii}, '_');
    roinames = names{3};

    
    atmp = splitroibeta(tmproibeta, tasknsub2);
%     ctmp = atmp;
%     atmp = {ctmp{1}, ctmp{3}, ctmp{2}, ctmp{5}, ctmp{4}};
    for mm = 1:5
        btmp{mm} = atmp{mm}/(std(atmp{mm})./sqrt(length(atmp{mm})));
        [~,p,~,stats] = ttest(btmp{mm});
        stmp.p(mm) = p;
        stmp.t(mm) = stats.tstat;
        stmp.df(mm) = stats.df;
    end
    roistats.(roinames) = stmp;
    
    meanforplot = cellfun(@mean, btmp);
    semforplot = cellfun(@(x) std(x)./sqrt(length(x)), btmp);
    mean_roi(ii,:) = meanforplot;
    sem_roi(ii,:) = semforplot;
    
    figure;
    taskcolors = [ 0 0 1;
        1 0 0; ... OC
        1 0 1; ... SC
        0 1 0; ... OA
        0,1,1]; % SA
    
    b = bar([1,2,3,4,5], meanforplot,'LineWidth', 1.2, 'BarWidth', 0.4);
    b.FaceColor = 'flat';
    b.CData = taskcolors;
    b.BaseLine.LineWidth = 1.3;
    hold on
    eb = errorbar([1,2,3,4,5], meanforplot, semforplot, 'o',...
        'LineWidth',2,'MarkerFaceColor','white','MarkerEdgeColor','black', 'CapSize',3, 'Color','black');
    box off
    
%     title(roinames)
    ylim([-4, 8])
    set(gca, 'ytick', -4:2:8, 'XTickLabel', {}, 'yticklabel', {}, 'linewidth', 1.3)

end