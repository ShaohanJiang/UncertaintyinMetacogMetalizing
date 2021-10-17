%% fig 2e
load('behtbl.mat')

nsub = size(behtbl, 1);
diffrt_crate_mtx = nan(4, 4, nsub);
% figure;
for jj = 1:nsub
    valididx = behtbl.resp{jj} >0 & behtbl.conf{jj} >0 ;
    thisisc = behtbl.isc{jj}(valididx);
    thisdiff = behtbl.diff{jj}(valididx);
    thisrt = behtbl.rtime{jj}(valididx);
    
    tmpprc = prctile( thisrt, [25 50 75]);
    tmpprcidx = thisrt < tmpprc;
    class1 = tmpprcidx(:,1);
    class2 = xor(tmpprcidx(:,1), tmpprcidx(:,2));
    class3 = xor(tmpprcidx(:,2), tmpprcidx(:,3));
    class4 = ~tmpprcidx(:,3);
    rtclassed = class1 + class2*2 + class3*3 + class4*4;
    
    for mm = 1:4
        for nn = 1:4
            thiscell = thisdiff ==mm & rtclassed==nn;
            thisrate = mean(thisisc(thiscell));
            diffrt_crate_mtx(mm, nn, jj) = thisrate;
            
        end
    end
    
%     subplot(4,7,jj)
%     heatmap(diffrt_crate_mtx(:,:,jj));
end

diffrt_crate_m = flipud(nanmean(diffrt_crate_mtx,3));
figure;
h = heatmap(diffrt_crate_m, 'Colormap', flipud(jet), 'ColorbarVisible',0);
h.ColorLimits = [0.25, 1];
h.XDisplayLabels= repmat({''},4,1);
h.YDisplayLabels= repmat({''},4,1);
% set(gca, 'YTickLabel', {}, 'XTickLabel', {})
% colorbar('off')
% cb = colorbar;
