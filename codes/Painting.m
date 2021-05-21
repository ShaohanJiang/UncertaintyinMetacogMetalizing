%% Painting

bncorr_bins = mean(after_corr, 3);
nsub = size(bncorr_bins,1);
mbins = mean(bncorr_bins);
sembins = std(bncorr_bins)/sqrt(nsub);
figure;
errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black')
% set(gca,'XTick',1:4, 'YTick',0.5:0.1:0.9, 'LineWidth',2.5, 'XTickLabel',{}, 'YTickLabel',{})
box off
xlim([0.5,10.5])
ylim([-0.1, 0.15])
hold on
plot([0.5, 10.5], [0, 0], 'k--','linewidth',1.2)
set(gca, 'XTick', 1:10, 'YTick',-0.05:0.05:0.15, 'LineWidth',1)

set(gca, 'XTickLabel',{}, 'YTickLabel',{});

atmp =mean(after_corr,3);
[h,p] = ttest(atmp(:,1),atmp(:,2)) % ttest n=2 and n=1 

%%
neauc_bins = mean(after_auc, 3);
nsub = size(neauc_bins,1);
mbins = mean(neauc_bins);
sembins = std(neauc_bins)/sqrt(nsub);
figure;
errorbar(mbins, sembins, 'o-','LineWidth',2.5,'MarkerFaceColor','white','MarkerEdgeColor','black')
% set(gca,'XTick',1:4, 'YTick',0.5:0.1:0.9, 'LineWidth',2.5, 'XTickLabel',{}, 'YTickLabel',{})
box off
xlim([0.5,10.5])
ylim([0.46, 0.54])
set(gca, 'XTick', 1:10, 'YTick',0.48:0.02:0.54, 'LineWidth',1);
hold on
plot([0.5, 10.5], [0.5, 0.5], 'k--','linewidth',1.2)
set(gca, 'XTickLabel',{}, 'YTickLabel',{});

%% paint confidence correlation before and after shuffle only
corr_shufflemean_bins = squeeze(mean(after_corruncer_resid))';
figure;
violinplot(corr_shufflemean_bins, [], 'ViolinColor', [0 0.4470 0.7410], 'ShowData', false, 'ViolinAlpha', 1, ...
    'EdgeColor', [0 0.4470 0.7410], 'BoxColor',  [0 0.4470 0.7410], 'ShowMean', true, 'MedianColor',[0 0.4470 0.7410]);

legend off
box off
xlim([0.5,10.5])
ylim([-0.5, 1])
set(gca, 'XTick', 1:10, 'YTick',-0.5:0.5:1, 'LineWidth',1);
hold on
plot([0.5,10.5], [0,0], 'k--')

hold on 
mbins = mean(corr_shufflemean_bins);
sembins = std(corr_shufflemean_bins); % range(CI)/2; %
errorbar(1:10, mbins(1:10), sembins(1:10), '.-','Color','k','LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});



%%
[nsub ,maxlen, ntimes ]= size(after_auc);

after_corrauc = zeros(ntimes, maxlen);

for jj = 1:ntimes
    after_corrauc(jj, :) = corr(before_auc, after_auc(:,:,jj));
    
end

%%
% figure; histogram(after_corrauc(:,10))
% box off
% hold on; plot([0,0], [0,700], '--', 'LineWidth',1.5)
figure;
violinplot(after_corrauc, [], 'ViolinColor', [0 0.4470 0.7410], 'ShowData', false, 'ViolinAlpha', 1, ...
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
mbins = mean(after_corrauc);
sembins = std(after_corrauc); %range(CI)/2; %
errorbar(1:10, mbins(1:10), sembins(1:10), '.-','Color','k','LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});

%%
meanauc = mean(after_auc,3);
figure;
violinplot(meanauc, [], 'ViolinColor', [0 0.4470 0.7410], 'ShowData', false, 'ViolinAlpha', 1, ...
    'EdgeColor', [0 0.4470 0.7410], 'BoxColor',  [0 0.4470 0.7410], 'ShowMean', true, 'MedianColor',[0 0.4470 0.7410]);
% violinplot(meanauc, [],  'ShowMean', true)
legend off
box off
set(gca, 'XTick',1:10)
xlim([0.5,10.5])
ylim([0.35, 0.7])
set(gca, 'XTick', 1:10, 'YTick',-0.4:0.1:0.7, 'LineWidth',1);
hold on
plot([0.5,10.5], [0.5,0.5], 'k--')

hold on 
mbins = mean(meanauc);
sembins = std(meanauc)/sqrt(size(meanauc,1)); %range(CI)/2; %
errorbar(1:10, mbins, sembins, '*','Color','k','LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black')
set(gca, 'XTickLabel',{}, 'YTickLabel',{});



%%
% after_corrauc = zeros(ntimes, maxlen);
% 
%     for jj = 1:ntimes
%         after_corrauc(jj, :) = corr(before_auc, after_auc(:,:,jj));
% 
%     end
% 
% figure; histogram(after_corrauc(:,10))
% box off
% hold on; plot([0,0], [0,700], '--', 'LineWidth',1.5)


%%
% figure;
% violin(after_corrauc, 'facecolor',[0 0.4470 0.7410],'edgecolor','none',...  % 1 0.56 0
% 'bw',0.3,'mc','k','medc',[], 'facealpha',0.75);
% legend off
% box off
% set(gca, 'XTick',1:10)
% xlim([0.5,10.5])
% ylim([-1.75, 0.54])
% set(gca, 'XTick', 1:10, 'YTick',0.48:0.01:0.54, 'LineWidth',2.5);
% hold on
% 
% plot([0,11], [0,0], 'k--')

%%
% p = [];
% CI = [];
% for ii = 1:size(after_corrauc,2)
%     tmp = randn(size(after_corrauc, 1),1)*std(after_corrauc(:,ii));
%     p(1,ii) = kstest2(after_corrauc(:,ii), tmp);
%     [~, p(2,ii), CI(:,ii)] = ztest(after_corrauc(:,ii), 0, std(after_corrauc(:,ii)));
%     
% end
