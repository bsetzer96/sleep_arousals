function [lagmean,lower,upper,cormean, lagstd] = bootstrapanal_heirarch_new(all_arousals, t, ROIs, totalArousAll)
%%
% load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_thal_20s.mat')
% thalnuc=allArousAll;
% thalrois=rois_hdr;
% load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_arous_20s.mat')
% thal=allArousAll(:,1,:);
% 
% all_arousals=thalnuc;
% all_arousals(:,10,:)=thal;
% ROIs=[thalrois, 'Thalamus'];

%%
% all_arousals=allArousAll;
% ROIs=rois_hdr;
%t=time;

%%
%ROIs=hdr.ROIs;
%t=hdr.t;
numrois=length(ROIs);

%dt=t(2)-t(1);
dt=0.247/4;
%range for which lags will be evaluated
range=50;

%cropping TS
ind= 60:(length(all_arousals)-60) ;%find((t>-10)&(t<10));
%preallocating space
Cor=zeros(numrois,numrois,1000);
Lcor=zeros(numrois,numrois,1000);
%e
thal_roi = all_arousals(:,1:numrois,:);
thal_ROIs = ROIs %{'PUL', 'MD', 'AV'};
subjArous=diff(totalArousAll);
%%
% n = # of arousals total
%n=112;
n=size(all_arousals,3);
% i = # of bootstraps
for i = 1:1000
    subind=ceil(rand(length(subjArous),1)*length(subjArous));
    for j=1:length(subind)
        k=subjArous(j);
        bind=ceil(rand(k,1)*k);
        boot_rois_subj=thal_roi(:,:,bind+totalArousAll(j));
        boot_rois_all(:,:,(totalArousAll(j)+1):totalArousAll(j+1))=boot_rois_subj;
    end
    boot_rois= mean(boot_rois_all,3); %averaging
    for c=1:numrois %for each ROI
        for d=1:numrois
            [xc, lags] = xcorr(boot_rois(ind,c),boot_rois(ind,d),range,'coeff');
            %figure(); plot(lags, xc); title([thal_ROIs(c) ' vs ' thal_ROIs(d)])
            [m,in]= max(abs(xc)); %extracting max correlation 
            Lcor(c,d,i)= lags(in)*dt; %defining lag of max cor
            Cor(c,d,i)= xc(in);
        end
    end
end
%%

lagmean = mean(Lcor,3);
%lagpct= prctile(
%lagstd = std(Lcor, [],3);


cormean=mean(abs(Cor),3);

bnds=prctile(Lcor, [2.5, 97.5],3);
lower=bnds(:,:,1);
upper=bnds(:,:,2);
lagstd=std(Lcor, [],3)
%%

avglag=nanmean(lagmean+diag(NaN*ones(length(ROIs),1)), 1)
%avglag=nanmean(lagmean, 1)
avgcor=nanmean(cormean+diag(NaN*ones(length(ROIs),1)), 1);
%%
up=nanmean(upper+diag(NaN*ones(length(ROIs),1)));
lw=nanmean(lower+diag(NaN*ones(length(ROIs),1)));


figure()
X=categorical(ROIs);
b = bar(X,avglag); 
xlabel('Regions'); ylabel('Lag (S)'); title('Average lags of averaged ROIs bootstrap')
hold on
er = errorbar(X,avglag,abs(avglag-lw), abs(up-avglag));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean Lag', '95% CI')
%% lag relative to whole thalamus
figure()
X=categorical(ROIs);
b = bar(X,lagmean(10,:)); 
xlabel('Regions'); ylabel('Lag (S)'); title('Lag of avg to thal bootstrap')
hold on
er = errorbar(X,lagmean(10,:),abs(lagmean(10,:)-lower(10,:)), abs(upper(10,:)-lagmean(10,:)));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean Lag', '95% CI')
%%
% n=size(all_arousals,3);
% Lcorj=zeros(numrois,numrois,n);
% Corj=Lcorj;
% for i = 1:1000
%     bind=ceil(rand(n,1)*n); %choosing random arousals
%     boot_rois_all = thal_roi(:,:,bind); %extracting random arousls
%     %boot_rois= mean(boot_rois_all,3); %averaging
%     for j=1:n %for each arousal
%         for c=1:numrois %for each ROI
%             for d=1:numrois
%                 [xc, lags] = xcorr(boot_rois_all(ind,c,n),boot_rois_all(ind,d,n),range,'coeff');
%                 %figure(); plot(lags, xc); title([thal_ROIs(c) ' vs ' thal_ROIs(d)])
%                 [m,in]= max(xc); %extracting max correlation 
%                 Lcorj(c,d,j)= lags(in)*dt; %defining lag of max cor
%                 Corj(c,d,j)= max(xc);
%             end
%         end
%     end
%         Lcor(:,:,i)=mean(Lcorj,3);
%         Cor(:,:,i)=mean(Corj,3);
% end

%%


end


