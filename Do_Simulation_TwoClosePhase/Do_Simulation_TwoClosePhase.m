%%
clear; home
fpara = {'rdm_osti_s1.txt', 'rdm_ostidmuncer_s1.txt', 'rdm_oconf_s1.txt', 'rdm_oconfdmuncer_s1.txt', ...
    'rdm_osti_s2.txt', 'rdm_ostidmuncer_s2.txt', 'rdm_oconf_s2.txt', 'rdm_oconfdmuncer_s2.txt'};
fsub = dir('s0*');

times = 1e3;

TR = 2;
nvolumn = 310; % 280 or 340
interval = 0.2;
hrf = FslHRF(0:interval: 36);

nsub = length(fsub);
npnt = 30;
nbeta = 10;
repeat = 100;
SNRrange = linspace(-30, 0, npnt);
% SNRrange = logspace(-40, 0, npnt);
% SNRrange = 10.^SNRrange;
% SNRrange = repmat(SNRrange, repeat, 1);
% SNRrange = SNRrange(:);

betahatb2 = zeros(4, nbeta, npnt, times, nsub);
b1all = randi([200,800], nbeta, npnt, nsub)/1000;
b3all = randi([200,800], nbeta, npnt, nsub)/1000;
b4 = 1*1e-5;  % set first phase effect to zero
b2all = randi([200,800], nbeta, npnt, nsub)/1000; % set second phase effect weight randomly

cor_twophase = zeros(nsub);

st = tic;
for jj = 1:nsub % for each sub
    
    % read behavior data 'rdm_osti_s1.txt, rdm_ostidmuncer_s1.txt,
    % rdm_oconf20_s1.txt, rdm_oconfdmuncer_s1.txt'
    ostimean = cell2mat( ReadText(fullfile(fsub(jj).name, fpara{1})));
    ostiuncer = cell2mat( ReadText(fullfile(fsub(jj).name, fpara{2})));
    oconfmean = cell2mat( ReadText(fullfile(fsub(jj).name, fpara{3})));
    oconfuncer = cell2mat( ReadText(fullfile(fsub(jj).name, fpara{4})));
    
    ostimat = ConvolveEV(ostimean, hrf, TR, nvolumn);
    ostiuncmat = ConvolveEV(ostiuncer, hrf, TR, nvolumn);
    oconfmat = ConvolveEV(oconfmean, hrf, TR, nvolumn);
    oconfuncmat = ConvolveEV(oconfuncer, hrf, TR, nvolumn);
    %%%[^^^] gen neural signals
    
    cor_twophase(jj)  = corr(ostiuncmat, oconfuncmat);
    
    for ss = 1:npnt % for each SNR
        thisSNR = SNRrange(ss);
        for mm = 1:nbeta
            b2 = b2all(mm, ss, jj);
            b1 = b1all(mm, ss, jj);
            b3 = b3all(mm, ss, jj);
            
            parfor tt = 1:times
                
                designmat = [ostimat, ostiuncmat, oconfmat, oconfuncmat];
                
                X = ostimat*b1 + ostiuncmat*b2 + oconfmat*b3 + oconfuncmat*b4;
                Y = awgn(X, thisSNR, 'measured');
                
                %         sigPower = sum(X.^2)/length(X)             %求出信号功率
                %         noisePower=sum((Y-X).^2)/length(Y-X)         %求出噪声功率
                %         SNR_10=10*log10(sigPower/noisePower)            %由信噪比定义求出信噪比，单位为db
                %         b=snr(X,Y-X)                             % snr(a,b) : a是原始信号，b是噪声信号
                %         figure;plot(1:length(X), [X, Y])
                
                designmat = bsxfun(@minus, designmat, mean(designmat));
                idmat    = pinv(designmat);
                betahat(:, mm, ss, tt, jj) = idmat * Y;
            end
        end
    end
end

% ttest(bhat(3,:), b1);
fprintf('DONE in %ss!\n',num2str(toc(st)))







