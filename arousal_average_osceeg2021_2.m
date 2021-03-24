

clc;clear;close all

cd /ad/eng/research/eng_research_lewislab/users/bsetzer/scripts/sleep_arousals/arousal_BOLD
addpath(genpath('../gen_functions'))

%version 2 is to extract ROIs out my my folder instead of labs data folder
%%
subjects = {'osceeg15b', 'osceeg17b', 'osceeg19b', 'osceeg20b', 'osceeg21b', 'osceeg22b'};
%subjects={'osceeg22b'};
%All cortical ROIs, except 'ctx-corpuscallosum'
 rois_hdr={'ctx-bankssts', 'ctx-caudalanteriorcingulate', 'ctx-caudalmiddlefrontal',...
    'ctx-cuneus','ctx-entorhinal', 'ctx-frontalpole', 'ctx-fusiform', 'ctx-inferiorparietal', 'ctx-inferiortemporal',...
    'ctx-insula', 'ctx-isthmuscingulate', 'ctx-lateraloccipital', 'ctx-lateralorbitofrontal', 'ctx-lingual',...
    'ctx-medialorbitofrontal', 'ctx-middletemporal', 'ctx-rostralmiddlefrontal', 'ctx-paracentral',...
    'ctx-parahippocampal', 'ctx-parsopercularis',  'ctx-parsorbitalis', 'ctx-parstriangularis', 'ctx-pericalcarine',...
    'ctx-postcentral', 'ctx-posteriorcingulate', 'ctx-precentral', 'ctx-precuneus', 'ctx-rostralanteriorcingulate', 'ctx-superiorfrontal',...
    'ctx-superiorparietal', 'ctx-superiortemporal', 'ctx-supramarginal', 'ctx-temporalpole', 'ctx-transversetemporal', 'ctx-unknown'};
%rois_hdr = {'wholeThalamus', 'brainstem', 'cortex'}
%rois_hdr = {'wholeThalamus',  'cortex'}
%rois_hdr = {'ctx_atlas', 'thalamus_atlas'};

%choose cortical folder for ctx and probinterp folder for thalamic nuclei
fold='cortical/';
%fold='probinterp/';
%fold='';

toplot=0;
%number of interpolated time points to use
range =217;%150; %300;%150;%300; %150;% 150;%300; %150; %300; % 150; %
%range in seconds for plots
rngsec=20;
%minimum and maximum number of clicks after arousal
minClick = 1 ;%1 ; %4
maxClick =25 ; %25

%matrix of how many arousals each subject has
totalArousAll = zeros(length(subjects),1);%%
%matrix of time x ROI x Arousal to store arousal timeseries
allArousAll = zeros(range*2+1, length(rois_hdr), totalArousAll(1));
figure()
for k = 1:length(subjects)
    subject = subjects{k};
    %datapath = ['/ad/eng/research/eng_research_lewislab/data/osceeg_frommgh/' subject];
    datapath = ['/projectnb/fastfmri/bsetzer/sleep_arousals/' subject];
    datapath2=['/ad/eng/research/eng_research_lewislab/data/osceeg_frommgh/' subject] ;
    fname=dir([datapath '/scores/' '*_score.mat']);
    path=[datapath '/rois/',fold];
    folder='osceeg_frommgh/';
    %folder='osceeg/';

    %storing number of arousals for each run
    totalArous = zeros(length(fname)+1,1);
    %storing arousal timeseries for each subject 
    allArous = zeros(range*2+1, length(rois_hdr), totalArous(1));
    for i = 1:length(fname) %loop through each run
        run = fname(i).name(1:5);
        cutoff = 20; %20; %number of seconds without a click to count as arousal
        %load info
        load([datapath2 '/eeg/fs200_averef_regbcg/' run '/hdr.mat'])
        timMRI=hdr.trs;
        %getting arousal times and clicks for run
        [arousalTimes, clicks]= extractArous(subject, run, cutoff, folder);
        %getting arousal timeseries and number of arousals
        [roiArous, time, kept, tracking] = extractArousalTS_osceeg2(subject, run, rois_hdr, folder, path, toplot, range, minClick, maxClick, arousalTimes, clicks, timMRI);
        %[roiArous, time, kept] = extractArousalTS(subject, run, rois_hdr, folder, path, toplot, range, minClick, maxClick, arousalTimes, clicks, timMRI);

        %store number of arousals
        %kept
        totalArous(i+1) = totalArous(i)+kept;
        if i ==1 
            %storing ROI timeseries matricies
            allArous(:,:,1:kept) = roiArous;
        else
            allArous(:, :, (totalArous(i)+1):totalArous(i+1)) = roiArous;
        end
    end

    %average for each subject accross runs
    allArousAvg = mean(allArous, 3);
    
    %figure();
    subplot(4,4,k)
    plot(time, allArousAvg); hold on
    plot([0 0], [min(allArousAvg, [], 'all') max(allArousAvg, [], 'all')])
    hold off
    xlabel('Time (S)')
    ylabel('% Change (BOLD)')
    title([subject ' Average ROIs at Arousal'])
    
    %stoing total number of arousals for subjects combined
    totalArousAll(k+1) = totalArousAll(k)+totalArous(end);
    if k ==1 
        allArousAll(:,:,1:totalArous(i+1)) = allArous;
    else
        allArousAll(:, :, (totalArousAll(k)+1):totalArousAll(k+1)) = allArous;
    end    
    
end
legend([ rois_hdr 'arousal'])
  %%
groupAvg = nanmean(allArousAll, 3);
groupStd= nanstd(allArousAll, [], 3)/sqrt(size(allArousAll,3));
%number of arousals for each subject
%totalArousAll(find(totalArousAll>0));

figure(6);    
%plot([0 0], [-2 5], 'k', 'LineWidth', 2); hold on  
%[3 8 12 19]
plot([0 0], [-1 1.2], 'g', 'LineWidth', 2); hold on
for l =   1:length(rois_hdr)
    subplot(3,3,l)
    xlim([-rngsec rngsec])
    nc=groupAvg(:,l)-mean(groupAvg(60:160,l)); %centered at baseline
    plot(time, nc, 'LineWidth', 2); 
    y2 = nc+groupStd(:,l);
    plot(time, y2, 'b')
    y1 = nc-groupStd(:,l);
    plot(time, y1, 'b')
    set(gca, 'FontSize',10)
    refline(0, 0);
    ylim([-2 1])
    title([ rois_hdr{l}])
    %title('Arousal-locked Activation')
end
legend('Arousal', rois_hdr{[1 2]}, 'standard error')
xlabel('Time (s)'); ylabel('BOLD signal (% change)')
%legend('arousal', 'Thalamus',  'Cortex', 'standard error', 'Location', 'Northeast') 
%legend('arousal', 'VLP', 'Anteroventral', 'standard error', 'Location', 'Northeast')
%legend('arousal', 'ROI signal', 'standard error', 'Location', 'NorthEast')
%% cross-correlation of averaged signals
[Lcor, Cor]= crosscoranal(groupAvg, time, rois_hdr, range);
lcorAvg= mean(Lcor);
%% bootstrap analysis of correlations and lags
%averages regions then takes cross-correlation
[lagmean, upper, lower,cormean, lagstd] = bootstrapanal(allArousAll, time, rois_hdr);

%% plotting mean lag and max correlation of bootsratp
figure(); heatmap(lagmean); set(gca, 'xdata', rois_hdr, 'ydata', rois_hdr); 
cmap=jet(601);
colormap(gca, cmap);
colorbar;
caxis([-2 2]);
figure()
heatmap(cormean); cmap=gray; colormap(gca,cmap)
set(gca, 'xdata', rois_hdr, 'ydata', rois_hdr);
%% average lag
avglag=mean(lagmean, 1);
avgcor=mean(cormean,1);
%save('ctx_lastavglags.mat', 'avglag', 'rois_hdr')
%save('ctx_avgcor.mat', 'avgcor', 'rois_hdr')
%% Plotting averages with error bars

%load('../../outputs/breathhold/breathhold_lags_thal.mat')
up=mean(upper);
lw=mean(lower);

figure()
X=categorical(rois_hdr);
b = bar(X,avglag); 
xlabel('Regions'); ylabel('Lag (S)'); title('Average lags of averaged ROIs bootstrap')
hold on
er = errorbar(X,avglag,abs(avglag-lw), abs(up-avglag));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean Lag', '95% CI')

%% avg correlations from bootstrap of arousals averaged and then cross-correlated
figure()
c=bar(X, avgcor)
xlabel('Regions of Interest'); ylabel('Cor (S)'); title('avereged ROIs bootstrap average max correlations')

%% Area-area lags 
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
xlabel('Regions'); ylabel('Lag (S)'); title('Lags')
hold on
er = errorbar(X,lm,abs(lm-lw), abs(up-lm));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean Lag', '95% CI')

%% Taking lags of individual arousals first and then averaging.
%singleArousalLags(allArousAll, time, rois_hdr)


%% saving

save('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/ctx_rois_20s_3T', 'allArousAll', 'rois_hdr', 'totalArousAll', 'groupAvg', 'groupStd')

%save('/projectnb/fastfmri/bsetzer/sleep_arousals/avg_ts/lags_thal_20s', 'lcorStd', 'lcorAvg')

