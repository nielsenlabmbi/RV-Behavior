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

message = 'Which eye was UNCOVERED durring experiment, Left/Dominant or Right/Weak? [l/r]: ';
eyeVisual = input(message, 's');


fnames=dir(['Z:\Ferret Behavior\RDK\lesion\testing postlesion\Covered eye\Left Covered\', animalId '\f*']);

data=[];

trialcount=1;
for f=1:length(fnames)

    load(fullfile('Z:\Ferret Behavior\RDK\lesion\testing postlesion\Covered eye\Left Covered\', animalId,fnames(f).name),'-mat')

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
StimLevels = round(rightData(:,5)*100); %making into intergers and rounding
StimLevels=unique(StimLevels)';

% Get unsided NumPos = number correct per coherence
for i=1:length(StimLevels)
    idx=find((round(rightData(:,5)*100)) ==StimLevels(i));
    NumPos(i) = sum(rightData(idx,4))';
end

%OutOfNum is the total trials run at each coherence value
trialsRun = ones(length(rightData(:,2)), 1);
for i=1:length(StimLevels)
    idx=find((round(rightData(:,5)*100)) == StimLevels(i));
   OutOfNum(i)=sum(trialsRun(idx))'; 
end



%%
N = 400; % 1000 is recomended but number of bootstraps to run to determine confidence intervals, can lower to 50 for fast checks

pDevBool = 'yes'; %optional, have as an extant varialbe to run goodness of fit simulations

PF = @PAL_CumulativeNormal;

addpath('C:\Rachel\Matlab\Palamedes1_11_9\Palamedes');

% %Get StimLevels (coherence leves)
% StimLevels = round(data(:,5)*100); %making into intergers and rounding
% StimLevels=unique(StimLevels)';
% 
% % Get unsided NumPos = number correct per coherence
% for i=1:length(StimLevels)
%     idx=find((round(data(:,5)*100)) ==StimLevels(i));
%     NumPos(i) = sum(data(idx,6))';
% end
% 
% %OutOfNum is the total trials run at each coherence value
% trialsRun = ones(length(data(:,6)), 1);
% for i=1:length(StimLevels)
%     idx=find((round(data(:,5)*100)) == StimLevels(i));
%    OutOfNum(i)=sum(trialsRun(idx))'; 
% end

%% Perfom fits for right sided stimuli
 [paramsValues,deltaCoherence, Dev, pDev, ci_thr, ci_s, ci_data, boot_data, t_boot, sl_boot] =  PsychCurve_sided_dots(StimLevels, NumPos, OutOfNum, animalId,N,pDevBool);

 %% save right side data
%
rightResults.files = fnames;
rightResults.data = data;
rightResults.rightData = rightData;
rightResults.ci_95s = boot_data;
rightResults.dCoh = deltaCoherence;
%rightResults.Thr75 = Threshold75;
%rightResults.Thr82 = Threshold82;
rightResults.pDev = pDev;
rightResults.Dev = Dev;
rightResults.eyeVisual = eyeVisual;
rightResults.numTrials = trialcount;
rightResults.numSessions = numSessions;
rightResults.PorCorrect = ProportionCorrectObserved;
rightResults.ci_thr = ci_thr;
rightResults.ci_s = ci_s;
rightResults.ci_data = ci_data;
rightResults.thr_boot = t_boot;
rightResults.slope_boot = sl_boot;
rightResults.paramsValues = paramsValues;

cd 'Z:\Ferret Behavior\RDK\lesion\testing postlesion\Covered eye\Left Covered\RightSide';
save(animalId, 'rightResults');
 
 
 %%  Repeat process for left side stimuli
 
 idxLeft = find(data(:,3)==-1);
leftData(:,1) = data(idxLeft,1);
leftData(:,2) =data(idxLeft,2);
leftData(:,3) =data(idxLeft,5);
leftData(:,4) =data(idxLeft,6);
leftData(:,5) =data(idxLeft,7);

StimLevels = round(leftData(:,5)*100); %making into intergers and rounding
StimLevels=unique(StimLevels)';

% Get unsided NumPos = number correct per coherence
for i=1:length(StimLevels)
    idx=find((round(leftData(:,5)*100)) ==StimLevels(i));
    NumPos(i) = sum(leftData(idx,4))';
end

%OutOfNum is the total trials run at each coherence value
trialsRun = ones(length(leftData(:,2)), 1);
for i=1:length(StimLevels)
    idx=find((round(leftData(:,5)*100)) == StimLevels(i));
   OutOfNum(i)=sum(trialsRun(idx))'; 
end

%% fit left side stimuli data
 [paramsValues,deltaCoherence, Dev, pDev, ci_thr, ci_s, ci_data, boot_data, t_boot, sl_boot] =  PsychCurve_sided_dots(StimLevels, NumPos, OutOfNum, animalId,N,pDevBool);

 %% save left side data
leftResults.files = fnames;
leftResults.data = data;
leftResults.leftData = leftData;
leftResults.ci_95s = boot_data;
leftResults.dCoh = deltaCoherence;
%leftResults.Thr75 = Threshold75;
%leftResults.Thr82 = Threshold82;
leftResults.pDev = pDev;
leftResults.Dev = Dev;
leftResults.eyeVisual = eyeVisual;
leftResults.numTrials = trialcount;
leftResultsnumSessions = numSessions;
leftResults.PorCorrect = ProportionCorrectObserved;
leftResults.ci_thr = ci_thr;
leftResults.ci_s = ci_s;
leftResults.ci_data = ci_data;
leftResults.thr_boot = t_boot;
leftResults.slope_boot = sl_boot;
leftResults.paramsValues = paramsValues;

cd 'Z:\Ferret Behavior\RDK\lesion\testing postlesion\Covered eye\Left Covered\LeftSide';
save(animalId, 'leftResults');

%%
ProportionCorrectObserved=NumPos./OutOfNum; 
StimLevelsFineGrain=(min(StimLevels):max(StimLevels)/1000:max(StimLevels));
ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
 
figure('name','Maximum Likelihood Psychometric Function Fitting');
plot(StimLevels,ProportionCorrectObserved,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',StimLevels);
axis([min(StimLevels) max(StimLevels) 0 1]);
hold on;
plot(StimLevelsFineGrain,ProportionCorrectModel,'g-','linewidth',4);
hold on; plot(StimLevels,.5*(ones(size(StimLevels))), 'r--')
errorbar(StimLevels, ProportionCorrectObserved, ((ci_data(2,:)) - (ci_data(1,:))), ((ci_data(3, :)) - (ci_data(2,:))))
xlabel('coherence level');
ylabel('proportion correct');
set(gca, 'fontsize',11);
title(animalId)
%ci_data is the 95% confidence intervals





% Saved for working with Erika's version of the code.
%This function expects a matrix called experiments. Where every row is a
%trial and first colum is the unsided stimulus and the second columns is if
%it was a correct answer.

%PsychCurve_general_ori(data(:,[6,5]), animalId)

%Sided data with outputs
 %[deltaAngle, pv, sl, pco, ci_thr, ci_s, ci_data] = PsychCurve_sidepref_ori(data(:,[6,5]) ,numSessions, animalId, 50)

%[StimLevel, NumPos, OutOfNum] = process_experiments_unsided(data(:,[6,5]));
 
% PsychCurve_sidepref_ori(experiments, numSessions,ferret,N)
% used in bootstrap to find confidence intervals


