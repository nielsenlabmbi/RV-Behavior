% Derived from dotStairs.m to capture PDS data for inputs into
% Palamedes psychometric curve fitting protocols. This version prepares
% data for sided, porportion to the right, analysis. Modifications to set up
% date from the lesion studies which have 4 stimuli and 2 answers.
% RV 10/3/2023


%% organize data 
clear all
close all
% up-date for each animal and eye visual 

message = 'What is the animalID? ';
animalId=input(message, 's');
%animalId='fbad6';

message = 'Which eye was UNCOVERED durring experiment, Left/Dominant, Right/Weak, None? [l/r/n]: ';
eyeVisual = input(message, 's');


fnames=dir(['Z:\Ferret Behavior\RDK\lesion\testing prelesion\', animalId '\f*']);

data=[];

trialcount=1;
for f=1:length(fnames)

    load(fullfile('Z:\Ferret Behavior\RDK\lesion\testing prelesion\', animalId,fnames(f).name),'-mat')

    %collect data from each trial
    
    if length(PDS.data) >20 %to avoid loading short trials, reduced for lesion study
    for i=1:length(PDS.data)
        data(trialcount,1)=f;
        %stair(i)=PDS.data{i}.stimulus.stair; %doesn't work right
        data(trialcount,2)=PDS.data{i}.stimulus.direction; %0 - side 1, 180 - side 2
        data(trialcount,3)=PDS.data{i}.stimulus.stimSide; % -1 = Left, 1 = Right
        
        if isfield(PDS.data{i}.stimulus,'dotCoherence')
            data(trialcount, 4)=1; % indicates cohherence is <1
            data(trialcount, 5)=PDS.data{i}.stimulus.dotCoherence;
        else
            data(trialcount, 4)=0; % indicates coherence is = 1, can be usd for sorting data later
            data(trialcount, 5)=1;
            %data(trialcount, 5)=PDS.conditions{i}.dotCoherence;
        end

        %correct
        if ~isfield(PDS.data{i}.pldaps,'goodtrial')
            data(trialcount, 6)=0;
        else
            data(trialcount, 6)=PDS.data{i}.pldaps.goodtrial;
        end
            %column 6, 0 for incorrect and 1 for correct trials

        trialcount=trialcount+1;
    end
    end
   end

numSessions = f;

% add sided coherence, and choice right for sided assesmend later
data(:,7)=data(:,5); %duplicate coherence column
data(:, 8)=data(:,6); %duplicate correctness colums
idx=(data(:,2)==0); % makes variable idx which makes a logical 1 or 0 from the directional data, 1 for dir 0 and 0 for dir 180
data(idx,7)=-1*data(idx,7); % makes coherence value negative (-x) for trials with 180 direction
data(idx, 8)=~data(idx,8); % applies the directional data in idx to the 1 or 0 correct trial data in 6 resulting did they go right


%% do simple stats

%sessionCoh = an array of the minimum coherence per session
session=unique(data(:,1))';

for i=1:length(session)
    idx=find((data(:,1)) ==session(i));
    sessionCoh(i) = min(data(idx,5)');
end
sessionCoh = sessionCoh';
mMinCoh = (sum(sessionCoh))/(length(sessionCoh));

% Now do one for percent correct per session
for i=1:length(session)
    idx=find((data(:,1)) ==session(i));
    perCorrect(i) = (sum(data(idx,6)))/(length(data(idx,6)));
end
perCorrect = perCorrect';
mPerCorrect = (sum(perCorrect))/(length(perCorrect));


% percent correct per session
%perCorS = an array of the percent correct per session
%want to seperate out data for stimSide (data(:,3)) and make matricies for
%responses by side (-1 left, 1 right), pull out direction of stim (:,2) and
%correctness (:,6) and coherence (:5)

%% generate dataset for stimuli displayed on the right
idxRight = find(data(:,3)==1);
rightData(:,1) = data(idxRight,1); %session 
rightData(:,2) = data(idxRight,2); %stimDirection
rightData(:,3) = data(idxRight,5); %coherence
rightData(:,4) = data(idxRight,6); % correct
rightData(:,5) =data(idxRight,7);% coherence side

%% Prep for fitting for Right side trials
%Get StimLevels (coherence leves)
rStimLevels = round(rightData(:,5)*100); %making into intergers and rounding
rStimLevels=unique(rStimLevels)';

% Get unsided NumPos = number correct per coherence
for i=1:length(rStimLevels)
    idx=find((round(rightData(:,5)*100)) ==rStimLevels(i));
    rNumPos(i) = sum(rightData(idx,4))';
end

%OutOfNum is the total trials run at each coherence value
trialsRun = ones(length(rightData(:,2)), 1);
for i=1:length(rStimLevels)
    idx=find((round(rightData(:,5)*100)) == rStimLevels(i));
   rOutOfNum(i)=sum(trialsRun(idx))'; 
end

% get percent correct per stim level for right side
 rPorCorrect = rNumPos./rOutOfNum;

 
 figure('name','Portion Correct Stim On Right Side');
plot(rStimLevels, rPorCorrect,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',rStimLevels);
axis([min(rStimLevels) max(rStimLevels) 0 1]);
ylabel('proportion correct');
hold on; plot(rStimLevels,.5*(ones(size(rStimLevels))), 'r--')
yyaxis right
plot (rStimLevels, rOutOfNum, 'g.','markersize',40);
xlabel('coherence level: - Left Direction, + Right Direction');
ylabel('number of trials');
set(gca, 'fontsize',11);
title(animalId, 'Stim On Right Side')
 
 %%  Repeat process for left side stimuli
 
 idxLeft = find(data(:,3)==-1);
leftData(:,1) = data(idxLeft,1);
leftData(:,2) =data(idxLeft,2);
leftData(:,3) =data(idxLeft,5);
leftData(:,4) =data(idxLeft,6);
leftData(:,5) =data(idxLeft,7);

lStimLevels = round(leftData(:,5)*100); %making into intergers and rounding
lStimLevels=unique(lStimLevels)';

% Get unsided NumPos = number correct per coherence
for i=1:length(lStimLevels)
    idx=find((round(leftData(:,5)*100)) ==lStimLevels(i));
    lNumPos(i) = sum(leftData(idx,4))';
end

%OutOfNum is the total trials run at each coherence value
trialsRun = ones(length(leftData(:,2)), 1);
for i=1:length(lStimLevels)
    idx=find((round(leftData(:,5)*100)) == lStimLevels(i));
   lOutOfNum(i)=sum(trialsRun(idx))'; 
end

 lPorCorrect = lNumPos./lOutOfNum;


figure('name','Portion Correct Stim On Left Side');
plot(lStimLevels, lPorCorrect,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',lStimLevels);
ylabel('proportion correct');
axis([min(lStimLevels) max(lStimLevels) 0 1]);
hold on; plot(rStimLevels,.5*(ones(size(rStimLevels))), 'r--')
yyaxis right
plot (lStimLevels, lOutOfNum, 'g.','markersize',40);
xlabel('coherence level: - Left Direction, + Right Direction');
ylabel('number of trials');
set(gca, 'fontsize',11);
title(animalId, 'Stim On Left Side')

%% Evolution of Sessions
figure
plot(data(:,5))
xlabel('Coherence');
ylabel('Trial Number');
title(animalId, 'General Evolution of Coherence W/Left Eye Covered')


%unsided data across all trials
cohLevels = round(data(:,5)*100); %making into intergers and rounding
cohLevels=unique(cohLevels); % for some reason "unique" gives duplicate rows for 0.05 to 0.20. try 
correctTrial = data(:,6);

for i=1:length(cohLevels)
    idx=find((round(data(:,5)*100)) == cohLevels(i));
   propCorr(i)=sum(correctTrial(idx))/length(idx); 
end


figure
plot(cohLevels, propCorr);
xlabel('Coherence');
ylabel('Portion Correct');
title(animalId, 'Unsided Total Portion Correct')

% 
%Sided analysis
cohLevelsSided=round((data(:,7)*100));
cohLevelsSided=unique(cohLevelsSided);
RsideTrial = data(:,8);

for i=1:length(cohLevelsSided)
    idx=find((round(data(:,7)*100)) ==cohLevelsSided(i));
    propCorrSide(i) = sum(RsideTrial(idx))/length(idx);
end

for i=1:length(cohLevelsSided)
    idx=find((round(data(:,7)*100)) == cohLevelsSided(i));
   trialNumSided(i)=sum(RsideTrial(idx)); 
end

propCorrSide = propCorrSide.'; %i dont know why it comes out going the wrong way but this fixes it


figure
%plots to both negative and positive sides
scatter(cohLevelsSided, (propCorrSide*100)')
xlabel('Coherence');
ylabel('Percent of Trials');
title(animalId, 'Trials to Side 1 by Direction')

figure
plot( cohLevelsSided, propCorrSide,'o')
xlabel('Coherence');
ylabel('Portion of Trials to Spout 1');
title(animalId, 'Portion Correct by Side')

