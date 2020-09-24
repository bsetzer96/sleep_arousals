function [roiArous, time, kept] = extractArousalTS(subject, run, rois_hdr, folder, path, toplot, range, minClick, maxClick, arousalTimes, clicks)
% new version of loading individual subjects, plotting, and averaging their arousals
%INPUT :: subject - string of subject name ex) 
%           run - string of run number
%           rois_hdr - array of strings corresponding to thalamus roi
%           folder - subject type folder ex)
%           path - path to data
%           toplot = 1 if you want plots
%OUTPUT :: roiArous - 3D matrix of each ROI TS for each arousal
%           kept - number of arousals
%% Extracting rois


%remember start time for mghthal is at FIRST TR, FIRST TR= TIME 0
%[arousalTimes, clicks]= extractArous_mghthal(subject, run, cutoff, folder);
%[arousalTimes, clicks]= extractArous_mghthal(subject, run, cutoff, folder);

%load ROIs
numROI= length(rois_hdr);
ROIs = load([path rois_hdr{1} '_' run '_timecourse.txt']);

time_MRI = 0:0.247:(0.247*(length(ROIs)-1));
if toplot==1
    figure(); plot(time_MRI, ROIs- mean(ROIs)); hold on
end
for i = 2:length(rois_hdr)
    ROIs = [ROIs, load([path rois_hdr{i} '_' run '_timecourse.txt']);];
    if toplot ==1 
        plot(time_MRI, ROIs(:,i) - mean(ROIs(:,i))); 
    end
end
if toplot ==1
plot(clicks, 40*ones(length(clicks), 1), '*');
legend([rois_hdr, 'clicks']); 
xlabel('Time'); ylabel('BOLD Signal')
axis([0 2000 -100 100])
title([ subject ' ' run])
hold off
end
%% load movement data

time_MRIfd=time_MRI(2:end);
fd = framewisedisplacement(subject, run, 0);
if toplot ==1
    figure(); plot(time_MRIfd, fd)
    title('Movement')
end

%% interpolate 
maxTime = length(time_MRI)*4;
timeInterp = linspace(time_MRI(1), time_MRI(end), maxTime);
roiInterp= zeros(maxTime, numROI);
for j=1:numROI
    roiInterp(:,j) = interp1(time_MRI', ROIs(:,j), timeInterp);
end

%% find time of closest TR to arousal
dt = timeInterp(2)-timeInterp(1);
%throw out 'arousals' that occured after MR stopped
arousalTimes = arousalTimes(find(arousalTimes<timeInterp(end)));
numArous= length(arousalTimes);
%matrix of timeseries in arousal window
roiArous = zeros(range*2+1, numROI, 0);
%clicks in arousal window
clickTimes= zeros(numArous,1);
int = (range)*dt;
time = -int:dt:int;
arous_sum = 0;
all_arousals = zeros(range*2+1, numROI, arous_sum);
kept = 0; thrown_out = 0;


%for each arousal, make sure we can take a range around the arousal (not
%too close to beginning or end), Make sure there's not too much movement
%around it. Otherwise thrown out.

for j = 1:numArous
    roiCrop = zeros(range*2+1, numROI, numArous);
    arousTime = arousalTimes(j);
    clickind = find((timeInterp < arousTime +(1/2)*dt) & (timeInterp > arousTime - (1/2)*dt));
    timeClick = timeInterp(clickind);
    if (clickind+range > length(roiInterp)) || (clickind-range<0) %if range around arousal is out of range of vector
        thrown_out = thrown_out+1;
    else %if inside the range
        roiCrop = roiInterp(clickind-range:clickind+range, :);
        timeCrop = timeInterp(:,clickind-range:clickind+range);
        %time of button presses within range
        clicksCrop = clicks(find((clicks<max(timeCrop))&(clicks>min(timeCrop))));
        %analyze movement around click
        fdInd = find((time_MRIfd>timeCrop(1)) & (time_MRIfd < timeCrop(end)));
        %throw out arousal if there is movement over 0.3
        if sum(find(fd(fdInd)>0.3))>0 %if too much movement
            thrown_out = thrown_out+1;
        elseif length(clicksCrop)< minClick %if there's not enough clicks
            thrown_out = thrown_out+1
        elseif length(clicksCrop)> maxClick % if theirs too many clicks
            thrown_out = thrown_out+1
        else
            %crop and save each roi
            for k = 1:numROI
                roiC = roiCrop(:,k);
                %demean
                roiCropPC(:,k) = ((roiC - mean(roiC))/mean(roiC))*100;
            end
            kept = kept+1;
            roiArous(:,:,kept) = roiCropPC;
            clickTimes(j) = timeClick;
            if toplot ==1
                figure()
                for l = 1:length(rois_hdr)
                    subplot(2,2,l)
                    plot(timeCrop, roiCropPC(:,l)); hold on
                    plot([timeClick timeClick], [min(roiCropPC, [], 'all') max(roiCropPC, [], 'all')])
                    plot(clicksCrop, 5*ones(length(clicksCrop)), '*'); hold off
                    legend([rois_hdr{l}, 'Arousal', 'Clicks'])
                    title([ subject ' ' run ' ' rois_hdr{l} ' at arousal'])
                end
            end
        
        end
    end

end




