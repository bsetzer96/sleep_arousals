%% comparing tSNR and other stats
%2/25/21
%load data
load('arousal_ts_data/thalnuc_lags.mat')
load('arousal_ts_data/thalnuc_stats.mat')
load('arousal_ts_data/thalnuc_tSNR.mat')

%% plot tSNR by amp

[c,b]=sort(a, 'descend');
atsnr=roiSNR(b);
figure(); bar(atsnr)
xticklabels(rois_hdr(b))
ylabel('tSNR')
xlabel('thalamic nuclei (highest amp->lowest amp)')
title('tSNR of thalamic nuclei')

%% amplitude lags
corrcoef(a,lcorAvg)

figure(); plot(a,lcorAvg, '*')
xlabel('Amplitude')
ylabel('Lags')

%% amplitude tSNR

corrcoef(a, roiSNR)

figure(); plot(a,roiSNR, '*')
xlabel('Amplitude')
ylabel('tSNR')

%% tSNR lags

corrcoef(lcorAvg, roiSNR)

figure(); plot(lcorAvg,roiSNR, '*')
xlabel('Lags')
ylabel('tSNR')