%% physio phase during arousal
%2/25/21
cd /ad/eng/research/eng_research_lewislab/users/bsetzer/scripts/sleep_arousals/arousal_BOLD
clc;clear;close all
addpath(genpath('../gen_functions'))

%%
% load subject
subjects={'mghthalamus1', 'mghthalamus2', 'mghthalamus3', 'mghthalamus5', 'mghthalamus6', 'mghthalamus9', 'mghthalamus10', 'mghthalamus12', 'mghthalamus13'};
range=20; %seconds
allRespArous=[];
allHbArous=[];
allHbPhase=[];

allHbAmp=[];
allRespPhase=[];
totalArousAll=zeros(length(subjects)+1,1);
for k = 1:length(subjects)
    subject = subjects{k};
    datapath = ['/projectnb/fastfmri/bsetzer/sleep_arousals/' subject '/'];
    fname=dir([datapath 'stcfsl_mc2_tvreg/' '*.nii']);
    path=[datapath 'physiophase/'];
    if subject(1)=='m'
        folder='mghthalamus/';
    else
        folder='inkrefcap/';
    end

    %storing number of arousals for each run
    totalArous = zeros(length(fname)+1,1);
    %storing arousal timeseries for each subject 
    allArous = zeros(range*2+1, 2, totalArous(1));
    kept=0;
    totalArous=zeros(length(fname));
    for i = 1:length(fname) %loop through each run
        run = fname(i).name(1:5);
        cutoff = 20; %20; %number of seconds without a click to count as arousal
        %getting arousal times and clicks for run
        fd = framewisedisplacement(subject, run, 0);
        timeMRIfd=0:0.247:(length(fd)-1);
        [arousalTimes, clicks]= extractArous_mghthal(subject, run, cutoff, folder);
        %getting physiotimeseries and number of arousals
        load([path run '_hbresp.mat'])
        hbf=bandpass(hbf, [0.7, 1.8], 1000);
        timePys=(1:length(phr))/1000; %in seconds
        %roiArous, time, kept] = extractArousalTS(subject, run, rois_hdr, folder, path, toplot, range, minClick, maxClick, arousalTimes, clicks);
        for j=1:length(arousalTimes)
            arousalTime=arousalTimes(j);
            physInd= find((timePys < arousalTime+range)&(timePys > arousalTime-range));
            fdInd=find((timeMRIfd < arousalTime+range)&(timeMRIfd > arousalTime-range));
            if (length(physInd)==40000) && (max(fdInd)<length(fd))
                if (sum(find(fd(fdInd)>0.3))==0)
                    kept=kept+1;
                    respArous=resp(physInd);
                    respHilb=hilbert(respArous);
                    hbfArous=hbf(physInd)-mean(hbf(physInd));
                    hbfHilb=hilbert(hbfArous);
                    hbfAmp=abs(hbfHilb);
                    %hbfArous1=max(hbfArous(10000:15000))-min(hbfArous(10000:15000));
                    %hbfArous2=max(hbfArous(25000:30000))- min(hbfArous(25000:30000));
                    phcpArous=phcp(physInd);
                    phrArous=phr(physInd);
                    %allHbAmp(:, totalArousAll(k)+kept)=[hbfArous1; hbfArous2];
                    allHbAmp(:, totalArousAll(k)+kept)=hbfAmp;
                    allHbPhase(:,totalArousAll(k)+kept)= phrArous;
                    allHbArous(:,totalArousAll(k)+kept)= hbfArous;
                    allRespPhase(:,totalArousAll(k)+kept)= phcpArous;    
                    allRespArous(:,totalArousAll(k)+kept)=respArous;
                end
            end
        end
    end
    totalArousAll(k+1)=totalArousAll(k)+kept;
end
%%
    timePysArous=(-20000:(20000-1))/1000;
    n=size(allRespArous,2);

    %average for each subject accross runs
    allRespAvg = mean(allRespArous,2);
    allRespStd=std(allRespArous,[],2)/sqrt(n);
    figure(); subplot(1,2,1); 
    plot([0 0], [-0.015 0.015], 'g'); hold on
    plot(timePysArous, allRespAvg); 
    plot(timePysArous, allRespAvg-allRespStd, 'b--')
    plot(timePysArous, allRespAvg+allRespStd, 'b--')
    title('Average Raw Respiration')
        
%     allRespPhaseAvg = mean(allRespPhase,2);
%         allRespPhaseErr = std(allRespPhase,[],2)/sqrt(n);
%     subplot(5,1,2); plot(timePysArous, allRespPhaseAvg); hold on
%     plot(timePysArous, allRespPhaseAvg+ allRespPhaseErr, 'b--')
%     plot(timePysArous, allRespPhaseAvg- allRespPhaseErr, 'b--')
%     title('Average Respiratory Phase')
%     
%     allHbAvg = mean(allHbArous,2);
%     allHbErr = std(allHbArous,[],2)/sqrt(n);
%     subplot(5,1,3); plot(timePysArous, allHbAvg); hold on
%     plot(timePysArous, allHbAvg+ allHbErr, 'b--')
%     plot(timePysArous, allHbAvg- allHbErr, 'b--')
%     title('Average Heartbeat Phase')
%     
%         allHbPhaseAvg = mean(allHbPhase,2);
%     allHbErr = std(allHbPhase,[],2)/sqrt(n);
%     subplot(5,1,4); plot(timePysArous, allHbPhaseAvg); hold on
%     plot(timePysArous, allHbPhaseAvg+ allHbErr, 'b--')
%     plot(timePysArous, allHbPhaseAvg- allHbErr, 'b--')
%     title('Average Heartbeat Phase')
    
    allHbAmpAvg = mean(allHbAmp,2);
    allHbAmpErr = std(allHbAmp,[],2)/sqrt(n);
    subplot(1,2,2); 
    plot([0 0], [0.7*10^-3 2.3*10^-3], 'g'); hold on
    plot(timePysArous, allHbAmpAvg); 
    plot(timePysArous, allHbAmpAvg-allHbAmpErr, 'b--')
    plot(timePysArous, allHbAmpAvg+allHbAmpErr, 'b--')
    title('Average Heartbeat Amplitude')
    ylim([0.7*10^-3 2.3*10^-3])
 
%     subplot(4,1,4); boxplot(allHbAmp');
%      ttest(allHbAmp(1,:), allHbAmp(2,:))
%      title('Heartbeat Amplitude Before and After Arousal')
%      xticklabels({'before', 'after'})
 
 %paired ttest for amplitude before and after not significant difference