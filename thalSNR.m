%function tSNR = thalSNR(tstart, time, V)
%used 2/25/21
cd /ad/eng/research/eng_research_lewislab/users/bsetzer/scripts/sleep_arousals
subjects = {'mghthalamus1' 'mghthalamus2' 'mghthalamus3' ...
    'mghthalamus5' 'mghthalamus6' 'mghthalamus9' 'mghthalamus10'...
    'inkrefcap5', 'inkrefcap7b', 'inkrefcap8', 'inkrefcap11'...
    'mghthalamus12', 'mghthalamus13'};
%'lp'
%rois_hdr = {'pulvinar', 'MD', 'av', 'lgn', 'cem' ,'cl', 'cm', 'l-sg', 'ld', 'mdl', 'mdm', 'mgn', 'mv(re)', 'pc', 'pf', 'pt', 'pua', 'pui', 'pul', 'pum', 'va', 'vamc', 'vla', 'vlp', 'vm', 'vpl'};
%rois_hdr = {'pulvinar', 'MD', 'av', 'lgn'};
rois_hdr = {'pulvinar', 'MD', 'av', 'lgn', 'cm',  'va', 'vla', 'vlp','vpl'};
% rois_hdr={'ctx-bankssts', 'ctx-caudalanteriorcingulate', 'ctx-caudalmiddlefrontal',...
%     'ctx-cuneus','ctx-entorhinal', 'ctx-frontalpole', 'ctx-fusiform', 'ctx-inferiorparietal', 'ctx-inferiortemporal',...
%     'ctx-insula', 'ctx-isthmuscingulate', 'ctx-lateraloccipital', 'ctx-lateralorbitofrontal', 'ctx-lingual',...
%     'ctx-medialorbitofrontal', 'ctx-middletemporal', 'ctx-ostralmiddlefrontal', 'ctx-paracentral',...
%     'ctx-parahippocampal', 'ctx-parsopercularis',  'ctx-parsorbitalis', 'ctx-parstriangularis', 'ctx-pericalcarine',...
%     'ctx-postcentral', 'ctx-posteriorcingulate', 'ctx-precentral', 'ctx-precuneus', 'ctx-rostralanteriorcingulate', 'ctx-superiorfrontal',...
%     'ctx-superiorparietal', 'ctx-superiortemporal', 'ctx-supramarginal', 'ctx-temporalpole', 'ctx-transversetemporal', 'ctx-unknown'}
%rois_hdr={'wholeThalamus', 'cortex', 'brainstem'};
%fold='';
fold='probinterp';
%fold='cortical/';

range = 150;
tr = 0.247
startTimes = [ 100 25 ...
    50 25 ...
    25 25 25 ...
    25 400 ...
    25 50 25 ...
    25 25 ...
    25 25 25 ...
    50 1350 ...
    100 850 ...
    50 1280 ...
    50 50 ...
    50 ...
    50 50];

a=0;
tSNRall = zeros(length(startTimes), length(rois_hdr));
allArousAll = zeros(range+1, length(rois_hdr), 17);
for k = 1:length(subjects)
    subject = subjects{k};
    datapath = ['/projectnb/fastfmri/bsetzer/sleep_arousals/' subject '/'];
    fname = dir([datapath 'stcfsl_mc2_tvreg/' '*.nii']);
    folder = 'mghthalamus/';
    path = [datapath 'rois/' fold];
    for i = 1:length(fname) %loop through each run
        run = fname(i).name(1:5);
        a=a+1;
        tstart = startTimes(a);
        for j = 1:length(rois_hdr)
            roiTS = load([path '/' rois_hdr{j} '_' run '_timecourse.txt']);
            time = 0:tr:(tr*length(roiTS));
            startInd = find((time <(tstart+tr)) & (time>(tstart-tr)));
            endInd = startInd+range;
            snrV= roiTS(startInd:endInd);
            tSNR = mean(snrV)./std(snrV);
            tSNRall(a,j) = tSNR;
            allArousAll(:,j,a) = snrV;
        end
    end
end

%%
roiSNR = mean(tSNRall,1)
roiSTE = std(tSNRall,[],1)/sqrt(a);
figure();
bar(1:length(rois_hdr),roiSNR)
%set(gca, 'XTickLabels', rois_hdr)
text(1:length(rois_hdr), roiSNR, rois_hdr, 'FontSize', 15)

xlabel('Thalamic Nuclei'); ylabel('tSNR'); title('tSNR of Thalamic Nuclei')
hold on
er = errorbar([1:length(rois_hdr)],roiSNR,roiSTE,roiSTE);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean tSNR', 'Standard Error')

%% Bootstrap analysis of lags

t = 0:tr:tr*(range);
maxTime = length(t)*4;
timeInterp = linspace(t(1), t(end), maxTime);
roiInterp= zeros(maxTime, length(rois_hdr), size(allArousAll,3));
for i = 1:size(allArousAll, 3)
    for j=1:length(rois_hdr)
        roiInterp(:,j) = interp1(t', allArousAll(:,j, i), timeInterp);
    end
end
%%
[lagmean, upper, lower]=bootstrapanal(roiInterp, timeInterp, rois_hdr)
%%
figure();
lags=abs(lagmean)
b = bar([1:3],lags);%{'MD->PUL', 'PUL-> AV', 'MD->AV'}, 
set(gca, 'xticklabels', {'MD->PUL', 'PUL-> AV', 'MD->AV'})
xlabel('Thalamic Nuclei'); ylabel('Lag (S)'); title('Lags Between Activated Thalamic Nuclei')
hold on
er = errorbar([1:3],lags,lower,upper);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
legend('Mean Lag', '95% CI')
set(gca, 'fontsize', 36)
