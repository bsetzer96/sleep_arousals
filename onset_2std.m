clc;clear
%load('arousal_ROIs.mat')
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_arous_20s.mat')
load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_arous_20s.mat')
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_20s_3T.mat');
% tr1=0.367;
% rng=218;
tr1=0.247;
rng=301;

baseind2=160;
%baseind2=107;
base_m = mean(groupAvg(60:160, :));
base_sd = std(groupAvg(60:160,:));

% base_m = mean(groupAvg(40:107, :));
% base_sd = std(groupAvg(40:107,:));
%%
% t = hdr.t;
% 
% %% smooth
% numROI = size(total_avg,2);
% %moving average
% smoothed_avg = zeros(size(total_avg));
% numsmooth = 20;
% %weights = [1/5 1/5 1/5 1/5 1/5]';
% for i = 1:numROI
%     for j = (numsmooth+1):length(smoothed_avg)-numsmooth
%         smoothed_avg(j,i) = sum(total_avg((j-numsmooth):(j+numsmooth),i)./ (numsmooth*2+1) );%.*weights);
%     end
% end
% 
% 
% total_avg = smoothed_avg;
% 
% 
% %%
% 
% m2 = 125;
% m1=20;
% bsln_ind = m1:m2;
% fprintf('Baseline will be from %f to %f seconds before arousal\n', [t(m1) t(m2)])
%n = size(total_avg, 2);
n= size(groupAvg,2);
%%
% base_m = mean(total_avg(bsln_ind, :));
% base_sd = std(total_avg(bsln_ind,:));


%%
%rois_hdr = {'pulvinar', 'MD', 'av' 'lgn', 'cm', 'va', 'vla', 'vlp', 'vpl'}
% ons = [];
% figure()
% t=-((length(groupAvg)-301)*0.247/4):(0.247/4):((length(groupAvg)-301)*0.247/4);
% for i = 1:9
%     roi = groupAvg(:,i);
%     mn = base_m(i);
%     sd = base_sd(i);   
%     on_ind = find((roi> (mn + 4*sd)) | (roi < (mn-4*sd)));
%     ind = on_ind>160;
%     ind = on_ind(ind);
%     onset = t(ind(1));
%     ons(i)=onset;
%     
%     subplot(3,3,i); 
%     plot(t, roi); hold on;
%     plot( [onset onset], [-2 2]); hold off;
%     title(rois_hdr{i})
%     
% end
%%
% [o, p]=sort(ons)
% rois_hdr{p}
%VPL MD VA PUL CM VLA VLP LGN AV
%PUL CM VLA VPL MD VA LGN VLP AV 2 std
%PUL VPL MD VA VLA VLP AV LGN 3 std
%MD VPL PUL CM VLA VA VLP AV LGN 4 std
%% amplitude
[a,b]=max(abs(groupAvg));
%a=[groupAvg(b(1),1), groupAvg(b(2),2)];
[~,c]=sort(b);
rois_hdr{c}

for i=1:n
    a(i)=groupAvg(b(i),i);
end

%VPL CM PUL MD VLP AV VA LGN VLA
%%
tr=tr1/4;
%max20=0.2*(a-base_m);
t=-((length(groupAvg)-rng)*tr):tr:((length(groupAvg)-rng)*tr);
% figure()
% ons20=[];
% ttp=[];
% for i = 1:length(rois_hdr)
%     lm=0.2*(groupAvg(b(i),1)-base_m(i));
%     mn=groupAvg(:,i);
%     if lm>0
%     onind=find(mn>lm);
%     elseif lm<0
%     onind=find(mn<lm);
%     end
%     onin=onind(find(onind>160));
%     onsetT=t(onin(1));
%     ons20(i)=onsetT;
%     
%     subplot(3,3,i); 
%     plot(t, mn); hold on;
%     plot( [onsetT onsetT], [-2 2]); hold off;
%     title(rois_hdr{i})
%     
%     %time to peak
%     ttp(i)=t(b(i))
%     
% end
%%
% [l, p]=sort(ons20)
% rois_hdr{p}
%CM VPL PUL VLP VLA LGN MD AV VA





%% full width half max

%half max 
%max50=0.5*(a-base_m);
% figure()
% fwhm=[];
% for i = 1:length(rois_hdr)
%     lm=0.5*(groupAvg(b(i),1)-base_m(i));
%     mn=groupAvg(:,i);
%     if lm>0
%     onind=find(mn>lm);
%     elseif lm<0
%     onind=find(mn<lm);
%     end
%     onin=onind(find(onind>160));
%     onT=t(onin(1));
%     offT=t(onin(end));
%     fwhm(i)=offT-onT;
%     
%     subplot(3,3,i); 
%     plot(t, mn); hold on;
%     plot( [onT onT], [-2 2]);
%     plot( [offT offT], [-2 2]); hold off;
%     title(rois_hdr{i})
%     
% end
%% onset time and ttp

tr=tr1/4;
%max20=0.2*(a-base_m);
t=-((length(groupAvg)-rng)*tr):tr:((length(groupAvg)-rng)*tr);
ons20=[];
ttp=[];
figure()
for i = 1:length(rois_hdr)
    lm=0.2*(a(i)-base_m(i));
    mn=groupAvg(:,i);
    if a(i)>0
    onind=find(mn>(base_m(i)+lm));
    else 
    onind=find(mn<base_m(i)+lm);
    end
    onin=onind(find(onind>baseind2));
    onsetT=t(onin(1));
    ons20(i)=onsetT;
    
    subplot(3,3,i); 
    plot(t, mn); hold on;
    plot( [onsetT onsetT], [-2 2]); hold off;
    title(rois_hdr{i})
    
    %time to peak
    ttp(i)=t(b(i))
    
end
%%
[l, p]=sort(ons20)
rois_hdr{p}
%CM VPL PUL VLP VLA LGN MD AV VA


%% full width half max

%half max 
%max50=0.5*(a-base_m);
figure()
fwhm=[];
for i = 1:length(rois_hdr)
    lm=0.5*(a(i)-base_m(i))
    mn=groupAvg(:,i);
    %if lm>0
    onind=find(mn>(base_m(i)+lm));
    %elseif lm<0
    %onind=find(mn<lm);
    %end
    onin=onind(find((onind>160) & (onind<450)));
    onT=t(onin(1));
    offT=t(onin(end));
    fwhm(i)=offT-onT;
    
    subplot(3,3,i); 
    plot(t, mn); hold on;
    plot( [onT onT], [-2 2]);
    plot( [offT offT], [-2 2]); hold off;
    title(rois_hdr{i})
    
end


