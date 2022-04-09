%%
load('beta_diff_rt_across_task.mat')
[ntask1, ntask2] = size(alldiffcbeta);
for jj =  1:ntask1
    for ii = 1:ntask2
        diffcbeta = alldiffcbeta{jj,ii};
        rtimebeta = allrtimebeta{jj,ii};
        if isempty(diffcbeta)
            continue
        end
        
        fig = figure;
        [r(jj,ii), p(jj,ii)] = PlotScatterLine(diffcbeta(:,1), diffcbeta(:,2), [0 0 0]);
%         set(gca,'XTick',0:0.2:1, 'YTick',0:0.2:1, 'LineWidth',1.2)
%         title(''); xlim([0, 1]); ylim([0, 1]);
%         set(gca, 'YTickLabel', {}, 'XTicklabel', {})
%         alpha(0.5)
        
        hold on
        [r2(jj,ii), p2(jj,ii)] = PlotScatterLine(rtimebeta(:,1), rtimebeta(:,2), [240,59,32]/255);
        hold on
        plot([0,0], [-0.2,1], 'Color',[0.5,0.5,0.5,0.2],'LineWidth',1.2)
        plot([-0.2,1], [0,0], 'Color',[0.5,0.5,0.5,0.2],'LineWidth',1.2)
        plot([-0.2,1], [-0.2,1], 'Color',[0.5,0.5,0.5,0.2],'LineWidth',1.2)
        set(gca,'XTick',-0.2:0.2:1, 'YTick',-0.2:0.2:1, 'LineWidth',1.2)
        title('');
        xlim([-0.2, 1]);
        ylim([-0.2, 1]);
        set(gca, 'YTickLabel', {}, 'XTicklabel', {})
        alpha(0.5)
    end
end