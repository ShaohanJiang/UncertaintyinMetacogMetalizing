%%
resfiles = dir('uncerbeta*');


for ii = 1:length(resfiles)
    atmp = load(fullfile(resfiles(ii).folder, resfiles(ii).name));
    atmp(1:3, :) = []; % remove voxel coordinates
    tmpmean = mean(atmp, 2);
    
    btmp = splitroibeta(atmp);
    ctmp = btmp;
    btmp = {ctmp{1}, ctmp{3}, ctmp{2}, ctmp{5}, ctmp{4}};

    %%% ttest
    ntask = length(btmp);
    thistest = ones(ntask);
    for jj = 1:ntask
        [~,p] = ttest(btmp{jj}, 0);
        thistest(jj,jj) = p;
        if jj == ntask
            break;
        end
        for mm = jj+1:ntask
            [~,p] = ttest2(btmp{jj}, btmp{mm});
            thistest(jj,mm) = p;
        end
    end
    allroitest{ii} = thistest;
    %%%%%%%%%%%%%
    meanforplot = cellfun(@mean, btmp);
    semforplot = cellfun(@std, btmp)./(cellfun(@length, btmp).^0.5);
    mean_roi(ii,:) = meanforplot;
    sem_roi(ii,:) = semforplot;
    
    figure;
    taskcolors = [
        0 0 1;
        1 0 0; ... OC
        1 0 1; ... SC
        0 1 0; ... OA
        0,1,1]; % SA
    
    b = bar([1,2,3,4,5], meanforplot,'LineWidth', 1.2, 'BarWidth', 0.4);
    b.FaceColor = 'flat';
    b.CData = taskcolors;
    b.BaseLine.LineWidth = 1.5;
    hold on
    eb = errorbar([1,2,3,4,5], meanforplot, semforplot, 'o',...
        'LineWidth',2,'MarkerFaceColor','white','MarkerEdgeColor','black', 'CapSize',3, 'Color','black');
    box off
    
end
