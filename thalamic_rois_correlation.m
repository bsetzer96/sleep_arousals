%calculating correlations between regions
addpath '/net/engnas/Research/eng_research_lewislab/users/bsetzer/scripts/sleep_arousals/arousal_BOLD'

%%


load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_thal_20s.mat')
thal=groupAvg;
thal_rois=rois_hdr;
allthal=allArousAll;

load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_arous_20s.mat')
avgs= groupAvg;
avg_rois=rois_hdr;
allavgthal=allArousAll(:,1,:);
thalavg=groupAvg(:,1);

allROIs= [thal, avgs];
all_hdr= [thal_rois, avg_rois(1)];

tr=0.247/4;
range=300;
time=-range*tr:tr:range*tr;

allArous=zeros(range*2+1, 10, 99);
allArous(:,1:9,:)=allthal;
allArous(:,10,:)=allavgthal;

[lagmean,lower,upper,cormean, lagstd] = bootstrapanal_heirarch_new(allArous, time, all_hdr, totalArousAll);
[Lcor, Cor]=crosscoranal([thal thalavg], time, all_hdr, range);

%% Rise

thalavg=avgs(:,2);
[c,x]=max(thalavg);


thalavg1=thalavg(100:x);
thal1=thal(100:x,:);
thalavg2=thalavg(x:end-100);
thal2=thal(x:end-100,:);

%% 1st half correlation
% [Lcor, Cor]=crosscoranal([thal, thalavg], time, all_hdr(1:end-1), range);
% 
% [Lcor1, Cor1]=crosscoranal([thal1, thalavg1], time(1:x), all_hdr(1:end-1), range);
% 
% [Lcor2, Cor2]=crosscoranal([thal2, thalavg2], time(x:end), all_hdr(1:end-1), range);
%%
Lcor=zeros(length(thal_rois),3);
Cor=zeros(length(thal_rois),3);
range=300;
for i=1:length(thal_rois)
    [xc, lags] = xcorr(thalavg,thal(:,i),range,'coeff');
    [m,in]= max(abs(xc));
    Lcor(i,1)= lags(in)*dt;
    Cor(i,1)= round(xc(range+1),2); 
    [xc, lags] = xcorr(thalavg1,thal1(:,i),range,'coeff');
    [m,in]= max(abs(xc));
    %figure(); plot( time(1:x), thalavg1, time(1:x), thal1(:,i))
    Lcor(i,2)= lags(in)*dt;
    Cor(i,2)= round(xc(range+1),2);     
    [xc, lags] = xcorr(thalavg2,thal2(:,i),range,'coeff');
    [m,in]= max(abs(xc));
    Lcor(i,3)= lags(in)*dt;
    Cor(i,3)= round(xc(range+1),2); 
    %figure(); plot( time(x:end), thalavg2, time(x:end), thal2(:,i))
end
%%



figure(); heatmap(Lcor', 'Colormap', jet); title('Lag Values Between ROIs'); 
%xlabel('To'); ylabel('From')
set(gca, 'yData', {'whole thalamus', 'whole thalamus upswing', 'whole thalamus down swing'}, 'xData', all_hdr(1:end-2))
%% heirarchical bootstrap

allArous1=allArous(1:x,:,:);
allArous2=allArous(x:end,:,:);

[lagmean,lower,upper,cormean, lagstd] = bootstrapanal_heirarch_new(allArous, time, all_hdr(1:end-1), totalArousAll)
[lagmean1,lower1,upper1,cormean1, lagstd1] = bootstrapanal_heirarch_new(allArous1, time, all_hdr(1:end-1), totalArousAll)
[lagmean2,lower2,upper2,cormean2, lagstd2] = bootstrapanal_heirarch_new(allArous2, time, all_hdr(1:end-1), totalArousAll)



%%

figure()
for i=1:length(thal_rois)
    subplot(length(thal_rois), 1, i)
    plot(time, thal(:,i)); hold on;
    plot(time, thalavg)
end
%%
figure(); heatmap(Cor(end-1:end,1:end-2), 'Colormap', jet); title('Lag Values Between ROIs'); 
xlabel('To'); ylabel('From')
set(gca, 'yData', all_hdr(end-1:end), 'xData', all_hdr(1:end-2))
xlabel('To'); ylabel('From')
set(gca, 'FontSize',16)
title('Correlation Between ROIs')
