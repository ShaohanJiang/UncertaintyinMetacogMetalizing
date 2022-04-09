function [r, p, sh, lh] = PlotScatterLine(x, y, c, titletxt)
% plot scatter fig for tow vector data with same length
% also plot the linear regression line
% [240,59,32]/255

if nargin<4
%     titletxt = 'Please enter the title!';
    titletxt = '';
end

[r,p] = corr(x,y);
txt = sprintf('%s\nr=%.2f, p=%.2f',titletxt, r,p);
% figure;
sh = scatter(x, y, 9, c,'filled');
hold on
pf=polyfit(x,y,1);
x1=-1.5:0.1:1.5;
y1=polyval(pf,x1);
lh = plot(x1,y1,'color',c, 'linewidth',1.2);
% xlim([0,1])
% ylim([0,1])
title(txt)

end