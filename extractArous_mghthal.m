function [arousalTimes bp] = extractArous_mghthal(subject, run, cutoff, folder)
%% Extract Click times


% if subject == 'inkrefcap6'
    % path = ['/ad/eng/research/eng_research_lewislab/users/bsetzer/data/inkrefcap6_run02_clicktimes.mat'];
% else
%     path = ['/ad/eng/research/eng_research_lewislab/data/', folder ,subject, '/behav/', run,'_clicktimes.mat'];
% end

% subject = 'inkrefcap8';
% run = 'run02';
% cutoff = 20;
% folder='inkrefcap/'

path = ['/ad/eng/research/eng_research_lewislab/users/bsetzer/data/', folder ,subject, '/behav/', run,'_clicktimes.mat'];

%% Find behavior arousal time

load(path)
 
%% Interclick intervals
%take difference between click rates y(n+1)- y(n)
if (subject(1)=='i') & (subject(end-1:end)~='11')
    bp=bp';
end

xn1 = [bp'; 0];
xn = [0; bp'];
clickinterval= xn1-xn;
clickinterval=clickinterval(2:end-1);

arous = find(clickinterval>cutoff);
arousalTimes = bp(arous+1); %time that corresponds to max.. bp(arous+1);
%sleepTimes=bp(arous);
end

