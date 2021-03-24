%% cross-covariance%
function [Lcor, Cor] = crosscoranal(total_avg, t, ROIs, range)
%load(aroual_ROIs)
%ROIs=hdr.ROIs
%t=hdr.t;

% ROIs=rois_hdr;
% t=time;
% total_avg=allArousAvg;

n=size(total_avg,2);
dt=t(2)-t(1);

ind= 60:(length(t)-60) ;%find((t>-10)&(t<10));
%ind=1:length(t);
Cor=zeros(n,n);
Lcor=zeros(n,n);

for c=1:n
    for d=1:n
        [xc, lags] = xcorr(total_avg(ind,c),total_avg(ind,d),range,'coeff');
         %figure();plot(lags*dt,xc, '*'); xlabel('Lag (S)'); ylabel('Correlation');
%         title([ROIs{c} ' vs ' ROIs{d}]); 
%         figure(); plot(t(ind), total_avg(ind,c), t(ind), total_avg(ind,d))
%         legend(ROIs{c}, ROIs{d});title('Raw Time-Series'); xlabel('Time (S)'); ylabel('BOLD Signal (%)')
        [m,in]= max(abs(xc));
        %lag with max correlation
        Lcor(c,d)= lags(in)*dt;
        %correlation at lag=0
        Cor(c,d)= round(xc(range+1),2); %round(xc(in),2);
    end
end
lags(range+1)

%%
figure(); heatmap(Lcor, 'Colormap', jet); title('Lag Values Between ROIs'); 
xlabel('To'); ylabel('From')
set(gca, 'xData', ROIs, 'yData', ROIs)
figure(); heatmap(Cor, 'Colormap', jet)
set(gca, 'xData', ROIs, 'yData', ROIs)
xlabel('To'); ylabel('From')
set(gca, 'FontSize',16)
title('Correlation Between ROIs')


figure(); bar(mean(Cor)); set(gca, 'xticklabel', ROIs)
title('Average Correlation of Regions'); xlabel('ROI'); ylabel('Correlation Strength')


%% t-test for significant difference of LGN
% thal_cor = Cor;
% lgn=4;
% reg=1:4;
% pvals=[];
% for j = [1 2 3]
%     v = reg(~ismember(reg,[j lgn]));
%     [h,p] = ttest2(thal_cor(v,lgn),thal_cor(v,j))
%     
%     pvals=[pvals p];
% end
% 
% pvals_corrected=pvals*3

%hypothesis test demonstrated that lgn not different from others

