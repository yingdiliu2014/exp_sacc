function [eyexy, successTrials, normalTrials] = separateTrials(subjectname)

% Input: subject name, so that we know which participant's data we are
% working with.

% Output: 
% 1. All eye position data (4 columns and one column of absolute time)
% 2. The starting sample and ending sample, and face stimuli number of all
% the 140 successful trials
% 3. Trials which do NOT contain blink samples during face presentation (eg
% due to temporary lost of track of the eye)
% 4. Ideal trials: Trials in which the data do not exceed an allowing area
% of 25 pixels.


%% 1. load relevant data from the datafiles
eyedatafile=[subjectname '_MSexp_eye.txt'];
eyedatatxt=fopen(eyedatafile);
eyedata=textscan(eyedatatxt, '%n%n%n%n%n%n%n%n%n%n%n%n%n'); % 13

bef_dur=(eyedata{1}); % before or after stimuli
stimnumber=(eyedata{3});
eyexy(:,1)=(eyedata{4}); % eyex1
eyexy(:,2)=(eyedata{5}); % eyey1
eyexy(:,3)=(eyedata{6}); % eyex2 
eyexy(:,4)=(eyedata{7}); % eyey2
timeloop=(eyedata{10});
% missing value: timetrial
eyexy(1,5)=timeloop(1);
for ii=2:length(timeloop)
    eyexy(ii,5)=eyexy(ii-1,5)+timeloop(ii);
end

datafile=[subjectname '_MSexp.txt'];
datatxt=fopen(datafile);
data=textscan(datatxt, '%n%n%n%n%n%n%n%n');

t_stimnumber=data{1};
t_sacc=data{2};


%% 2. detecting blink samples (or samples where eye tracker lost track)
blinkthreshold = [1920 1080];
blink = zeros(length(stimnumber),1);
for sample=1:length(stimnumber)
    if abs(eyexy(sample,1))>blinkthreshold(1) || abs(eyexy(sample,2))>blinkthreshold(2) || abs(eyexy(sample,3))>blinkthreshold(1) || abs(eyexy(sample,4))>blinkthreshold(2)
        blink(sample)=1;
    end;
end


%% 3. separating individual trials from the many eye data samples
totaltrialnumber=length(t_sacc); % how many trials altogether has this subject did
successtrialnumber=length(find(t_sacc==0)); % should be 140
startOfTrials=zeros(totaltrialnumber,1);
startOfTrials(1)=1; % the first sample is always the start of the first trial
trialsCounted=1;
for ii=2:length(stimnumber)
    if stimnumber(ii)~=stimnumber(ii-1)
        trialsCounted=trialsCounted+1;
        startOfTrials(trialsCounted)=ii;
    end
end
endOfTrials=[startOfTrials(2:end)-1; length(stimnumber)];
trialSamples=[startOfTrials endOfTrials];


%% 4. select successful trials. Identify sample number and face number for each of them. 
successTrialSamples=trialSamples(t_sacc==0,:);
successTrialFaces=t_stimnumber(t_sacc==0);
successTrialSamplesAll=[]; % put all the samples together
for ii=1:length(successTrialSamples)
    successTrialSamplesAll=[successTrialSamplesAll successTrialSamples(ii,1):successTrialSamples(ii,2)];
end
successTrials=[successTrialSamples successTrialFaces];

%% 5. abort abnormal trials
blinksamples=find(blink==1);
wierdo=[];
for ii=1:length(blinksamples)
    if ismember(blinksamples(ii),successTrialSamplesAll)
        wierdo=[wierdo blinksamples(ii)];
    end
end
% if the wierdo is in the pre-stimulus period then it's fine, if not,
% then it's abnormal. 
abnormal=intersect(find(bef_dur==1),find(blink==1));
record=[];
for ii=1:length(abnormal)
    sample=abnormal(ii);
    for jj=1:length(successTrialSamples)
        if ismember(sample,successTrialSamples(jj,1):successTrialSamples(jj,2))
            break
        end
    end
    record(ii)=jj;
end
abnormalTrials=unique(record);
normalTrials=setdiff(1:successtrialnumber, abnormalTrials);


%% 6. select ideal trials (optional)
%% (in which both eyes cleanly stay within a certain (seuildrift) range. 
% outrageousTrials=[]; 
% for tr=1:length(normalTrials)
%     theTrial=normalTrials(tr);
%     samplesOfThisTrial=[successTrialSamples(theTrial,1):successTrialSamples(theTrial,2)]';
%     trialXys=eyexy(samplesOfThisTrial,:);
%     % trials that contain out-of-bound samples 
%     for ii=1:size(samplesOfThisTrial)
%         % left eye distance of this sample 
%         diffx1=trialXys(ii,1)-wRect(1)/2; 
%         diffy1=trialXys(ii,2)-wRect(2)/2;
%         distdrift1=sqrt(diffx1*diffx1+diffy1*diffy1);
%         % right eye distance of this sample 
%         diffx2=trialXys(ii,3)-wRect(1)/2;
%         diffy2=trialXys(ii,4)-wRect(2)/2;
%         distdrift2=sqrt(diffx2*diffx2+diffy2*diffy2);
%         distdrift=max(distdrift1,distdrift2); 
%         if distdrift>seuildrift
%             outrageousTrials=[outrageousTrials, tr];
%             break
%         end
%     end
% end
% idealTrials=setdiff(normalTrials, outrageousTrials);

end