function [lagmean,lower,upper,cormean] = bootstrapanal(all_arousals, t, ROIs)
%all_arousals=allArousAll;
%ROIs=rois_hdr;
%t=time;
%%
%load('arousal_ROIs_all.mat')
%%
%ROIs=hdr.ROIs;
%t=hdr.t;
numrois=length(ROIs);

dt=t(2)-t(1);
%range for which lags will be evaluated
range=50;

%cropping TS
ind= 60:(length(t)-60) ;%find((t>-10)&(t<10));
%preallocating space
Cor=zeros(numrois,numrois,1000);
Lcor=zeros(numrois,numrois,1000);
%excluding LGN
thal_roi = all_arousals(:,1:numrois,:);
thal_ROIs = ROIs %{'PUL', 'MD', 'AV'};
%%
% n = # of arousals total
%n=112;
n=size(all_arousals,3);
% i = # of bootstraps
for i = 1:1000
    bind=ceil(rand(n,1)*n); %choosing random arousals
    boot_rois_all = thal_roi(:,:,bind); %extracting random arousls
    boot_rois= mean(boot_rois_all,3); %averaging
    for c=1:numrois %for each ROI
        for d=1:numrois
            [xc, lags] = xcorr(boot_rois(ind,c),boot_rois(ind,d),range,'coeff');
            %figure(); plot(lags, xc); title([thal_ROIs(c) ' vs ' thal_ROIs(d)])
            [m,in]= max(xc); %extracting max correlation 
            Lcor(c,d,i)= lags(in)*dt; %defining lag of max cor
            Cor(c,d,i)= max(xc);
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


end


