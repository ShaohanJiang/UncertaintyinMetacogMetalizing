function PlotErrorbar(tbl, doplotsigle)
%   ʵ�ֹ��ܣ�plot ����򲻷�������ݵ�barͼ��errorbar��ͼ�����ҿ��ԽϷ����
%   ����һЩͼ������
%   ���ݽṹ��֧�� table cell mat, �����ԭʼֵ��ÿһ��Ϊһ���������������Զ�����mean��sem
%   Argument pairs:
%
%       (c) William Jiang, 2019
%
if nargin<2
    doplotsigle = true;
end

if iscell(tbl)
    % nsub = size(mat,1);
    % matmean = mean(mat);
    % matsem = std(mat)/sqrt(nsub);
%     gcolor = [0 0 1; 1 0 0; 1 0 1; 0 1 0; 0 1 1];
    matmean = cellfun(@mean, tbl);
    gcolor = lines(length(matmean));
    nmean = size(matmean,1);
    
    matsem = bsxfun(@rdivide, cellfun(@std, tbl), sqrt(cellfun(@length, tbl)));
    % rowname = tbl.Row;
    
    if nmean == 1
%         figure;
        b= bar(matmean);
        b.FaceColor = 'flat';
        b.CData = gcolor;
        %     set(gca,'xticklabel', rowname);
        if doplotsigle
            hold on
            for ii = 1:length(tbl)
                xcoor = repmat(b.XData(ii), length(tbl{ii}),1);
                scatter(xcoor, tbl{ii}, '.', 'k')
            end
        end
        hold on
        % xlabel(rowname')
        errorbar(matmean, matsem, 'linestyle', 'none','LineWidth',1.5,'Color','k', ...
            'CapSize',5); %'MarkerEdgeColor','black','Marker','o', 'MarkerFaceColor','white',
        hold off
    else
        
%         figure;
        b= bar(matmean);
        
        %     set(gca,'xticklabel', rowname);
        
        hold on
        % xlabel(rowname')
        errorbar(matmean, matsem, 'linestyle', 'none'); hold off
        
        
    end
else
    datmean = mean(tbl);
    sem = std(tbl)/sqrt(size(tbl,1));
%     figure;
    bar(datmean)
    hold on;
    errorbar(datmean, sem, 'ko');
    
end
hold off;
end % END OF FILE