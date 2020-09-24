% framewise displacement

% for a given subject/run load motion logs and calculate framewise
% displacement in mm

function fd=framewisedisplacement(subj,rnum,ploton);


%if subj(1) == 'i'
 %       chanset = 'stcfsl_mc2';
  %      datapath=['/ad/eng/research/eng_research_lewislab/data/inkrefcap/' subj '/' chanset '/'];
% elseif subj(1) == 'o'
%     chanset='stcfsl_mc2_zpad';
%     datapath=['/ad/eng/research/eng_research_lewislab/data/osceeg_frommgh/' subj '/' chanset '/'];
% elseif subj(1) == 'm'
    chanset = 'stcfsl_mc2';
    datapath=['/projectnb/fastfmri/bsetzer/sleep_arousals/' subj '/' chanset '/'];
% end
fname=dir([datapath 'mclogs/' rnum '*_mc2.mcdat']);

a=load([datapath 'mclogs/' fname(1).name]);

mm=a(:,5:7);
rot=a(:,2:4);
rotmm=deg2rad(rot)*50; % power 2012

dd=diff([rotmm mm]);
fd=sum(abs(dd),2);

if ploton
figure
hold on
subplot(2,1,1);
plot(a(:,2:7));
subplot(2,1,2);
plot(sum(abs(fd),2));
title([' subject number: ' subj ]) ; 
end
