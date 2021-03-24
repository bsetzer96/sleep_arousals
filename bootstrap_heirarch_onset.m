%bootstrap 20% rise time
clc;clear
%load arousals 
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/all_thal_arous_20s')
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_arous_20s')
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_20s_3T')
%load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/thal_ctx_NREM_3T')
load('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/ctx_rois_20s_3T.mat')

% totalArousAll=nREMArousAll;
% allArousAll=arousNREMall;


%tr1=0.247;
tr1=0.367;
%rng=300;
rng=218;
%base1=60;
%base2=160;
 base1=40;
 base2=107;
%%
subjArous=diff(totalArousAll);
thal_roi=allArousAll;
tr=tr1/4;
t=-tr*rng:tr:tr*rng;
nsubj=length(subjArous);
n=size(allArousAll,3);
numROIs=length(rois_hdr);
allOnset=zeros(1000, numROIs);
boot_rois_all=zeros(size(thal_roi));
%% i = # of bootstraps
for i = 1:1000 %1000
   subind=ceil(rand(length(subjArous),1)*length(subjArous));
    for j=1:length(subind)
        k=subjArous(j);
        bind=ceil(rand(k,1)*k);
        boot_rois_subj=thal_roi(:,:,bind+totalArousAll(j));
        boot_rois_all(:,:,(totalArousAll(j)+1):totalArousAll(j+1))=boot_rois_subj;
    end
    boot_rois= mean(boot_rois_all,3); %averaging
    %calculate basilines
    base_m = mean(boot_rois(base1:base2, :));
    %calculate max
    %[a,b]=max(abs(groupAvg));
    [a,b]=max(abs(boot_rois));
    %a=[groupAvg(b(1),1), groupAvg(b(2),2)];
    %calculate 20% rise time
    
    %figure()
    ons20=zeros(1, numROIs);
    
    %figure()
    for l = 1:numROIs
            a(l)=boot_rois(b(l),l);
            max20(l)=0.2*(a(l)-base_m(l));
            lm=max20(l);
            mn=boot_rois(:,l);
           if a(l)>=0
            onind=find(mn>(base_m(l)+lm));
           else 
             onind=find(mn<(base_m(l)+lm));
           end
            onin=onind(find(onind>base2));
            if isempty(onin)
                onsetT=NaN;
            else
            onsetT=t(onin(1));
            end
            ons20(l)=onsetT;
% 
%             subplot(3,3,i); 
%             plot(t, mn); hold on;
%             plot( [onsetT onsetT], [-2 2]); hold off;
%             title(rois_hdr{i})
        
%             subplot(3,3,j); 
%             plot(t, mn); hold on;
%             plot( [onsetT onsetT], [-2 2]); hold off;
%             title(rois_hdr{j})
    end
    allOnset(i,:)=ons20;
end

%%
avgOnset=mean(allOnset);
bnds=prctile(allOnset, [2.5, 97.5],1);
lower=bnds(1,:);
upper=bnds(2,:);

[l, p]=sort(avgOnset)
rois_hdr{p}
%CM VPL LGN VLP PUL VLA MD VA AV
%VLP = 0.55 : -5.2

%% sort from earliest to latest
% [l, p]=sort(ons20)
% rois_hdr{p}
%CM VPL PUL VLP VLA LGN MD AV VA
