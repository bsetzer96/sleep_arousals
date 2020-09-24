function [lagmean,lower,upper,cormean] = singleArousalLags(all_arousals, t, ROIs)
% all_arousals=allArousAll;
% ROIs=rois_hdr;
% t=time;
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


n=size(all_arousals,3);
Lcor=zeros(numrois, numrois, n);
Cor=Lcor;
for j=1:n %for each arousal
    for c=1:numrois %for each ROI
        for d=1:numrois
            [xc, lags] = xcorr(thal_roi(ind,c,n),thal_roi(ind,d,n),range,'coeff');
            %figure(); plot(lags, xc); title([thal_ROIs(c) ' vs ' thal_ROIs(d)])
            [m,in]= max(xc); %extracting max correlation 
            Lcor(c,d,j)= lags(in)*dt; %defining lag of max cor
            Cor(c,d,j)= max(xc);
        end
    end
end
lagmean=mean(Lcor,3);
cormean=mean(Cor,3);

%%
figure(); heatmap(lagmean); set(gca, 'xdata', ROIs, 'ydata', ROIs); 
cmap=jet(601);
colormap(gca, cmap); title('Individual Arousal Lag Average')
colorbar;
caxis([-2 2]);
figure()
heatmap(cormean); cmap=gray; colormap(gca,cmap)
set(gca, 'xdata', ROIs, 'ydata', ROIs);
title('Individual Aroousal Correlation Average')
%%
bnds=prctile(Lcor, [2.5, 97.5],3);
lower=bnds(:,:,1);
upper=bnds(:,:,2);

lower=std(Lcor,[],3)/sqrt(n);
upper=std(Lcor,[],3)/sqrt(n);
%%
m={'pul-pul', 'pul-md', 'pul-av', 'pul-lgn', 'pul-cm', 'pul-va', 'pul-vla','pul-vlp', 'pul-vpl',...
    'md-pul', 'md-md', 'md-av', 'md-lgn', 'md-cm', 'md-va', 'md-vla', 'md-vlp','md-vpl',...
    'av-pul', 'av-md', 'av-av', 'av-lgn', 'av-cm', 'av-va', 'av-vla', 'av-vlp','av-vpl',...
    'lgn-pul', 'lgn-md', 'lgn-av', 'lgn-lgn', 'lgn-cm', 'lgn-va', 'lgn-vla', 'lgn-vlp','lgn-vpl',...
    'cm-pul', 'cm-md', 'cm-av', 'cm-lgn', 'cm-cm', 'cm-va', 'cm-vla', 'cm-vlp','cm-vpl',...
    'va-pul', 'va-md', 'va-av', 'va-lgn', 'va-cm', 'va-va', 'va-vla', 'va-vlp','va-vpl',...
    'vla-pul', 'vla-md', 'vla-av', 'vla-lgn', 'vla-vm', 'vla-va', 'vla-vla', 'vla-vlp','vla-vpl',...
    'vlp-pul', 'vlp-md', 'vlp-av', 'vlp-lgn', 'vlp-vm', 'vlp-va', 'vlp-vla', 'vlp-vlp', 'vlp-vpl',...
    'vpl-pul', 'vpl-md', 'vpl-av', 'vpl-lgn', 'vpl-vm', 'vpl-va', 'vpl-vla', 'vpl-vlp','vpl-vpl'};
lm=reshape(lagmean, [81,1]);
lw=reshape(lower, [81,1]);
up=reshape(upper, [81,1]);
figure()
X=categorical(m);
b = bar(X,lm); 
xlabel('Regions'); ylabel('Lag (S)'); title('Individual Arousal Lags Averaged')
hold on
er = errorbar(X,lm,lw, up);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean Lag', '95% CI')

end