%% PCA func

%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_thal_20s.mat')
load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_thal_60s.mat')

thal=groupAvg;
thal_rois=rois_hdr;
%%
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_ctx_20s.mat')
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/ctx_rois_20s_3T.mat')
load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_ctx_60s.mat')

ctx= groupAvg;
ctx_rois=rois_hdr;
%%

X= [thal, ctx];
X_label= [thal_rois, ctx_rois];

% X= [thal];
% X_label= [thal_rois];

% X= [ctx];
% X_label= [ctx_rois];

%7T
% tr=0.247/4;
% time=-300*tr:tr:300*tr;

tr=0.247/4;
time=-300*3*tr:tr:300*tr;


%3T
% tr=0.367/4;
% time=-217*tr:tr:tr*217


%%
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(X, 'NumComponents', 3);
EXPLAINED
%%
figure(); plot(time, SCORE(:,1)); hold on; plot(time, SCORE(:,2)) ;plot([0,0], [-6, 4]);  legend({'1', '2',  'arousal'}); xlabel('Sample'); ylabel('BOLD signal')
%plot(SCORE(:,3))
title('PCA of all areas')
%%
figure(); heatmap(COEFF(:,1:2), 'colormap', jet); set(gca, 'ydata', X_label)
title('Principal Component Coefficients')
%%
%addpath(genpath('/ad/eng/research/eng_research_lewislab/users/bsetzer/sleep_arousals/FastICA_25'))
%note: ICA looks like it just pulls out noise
[icasig, A, W] = fastica (X', 'numOfIC', 3);
%%
figure(); %plot(icasig([4 5 6], :)'); legend({'1', '2', '3', '4', '5'}); xlabel('Sample'); ylabel('BOLD signal')
plot(icasig')
figure(); heatmap(A); set(gca, 'ydata', X_label)
