% Derived from dotStairs.m to capture PDS data for inputs into
% Palamedes psychometric curve fitting protocols. This version preserves
% workflow before breaking into sided and unsided data organization.

animalId='fbae2';

fnames=dir(['Z:\Ferret Behavior\RDK\Amblyopia\stairs dominant\' animalId '\f*']);

data=[];

trialcount=1;
for f=1:length(fnames)

    load(fullfile('Z:\Ferret Behavior\RDK\Amblyopia\stairs dominant\',animalId,fnames(f).name),'-mat')

    %collect data from each trial
    
    if length(PDS.data) >64 %to avoid loading short trials, want 65 or greater
    for i=1:length(PDS.data)
        data(trialcount,1)=f;
        %stair(i)=PDS.data{i}.stimulus.stair; %doesn't work right
        data(trialcount,2)=PDS.data{i}.stimulus.direction; %0 - side 1, 180 - side 2
        
        if isfield(PDS.data{i}.stimulus,'dotCoherence')
            data(trialcount,3)=1;
            data(trialcount,4)=PDS.data{i}.stimulus.dotCoherence;
        else
            data(trialcount,3)=0;
            data(trialcount,4)=PDS.conditions{i}.dotCoherence;
        end

        %correct
        if ~isfield(PDS.data{i}.pldaps,'goodtrial')
            data(trialcount,5)=0;
        else
            data(trialcount,5)=PDS.data{i}.pldaps.goodtrial;
        end
            %column 5, 0 for incorrect and 1 for correct trials

        trialcount=trialcount+1;
    end
    end
end

%add sided coherence, and choice right
data(:,6)=data(:,4); %duplicate coherence column
data(:,7)=data(:,5); %duplicate correctness colums
idx=(data(:,2)==0); % makes variable idx which makes a logical 1 or 0 from the directional data, 1 for dir 0 and 0 for dir 180
data(idx,6)=-1*data(idx,6); % makes coherence value negative (-x) for trials with 180 direction
data(idx,7)=~data(idx,7); % applies the directional data in idx to the 1 or 0 correct trial data in 5 resulting did they go right

numSessions = f;
%This function expects a matrix called experiments. Where every row is a
%trial and first colum is the unsided stimulus and the second columns is if
%it was a correct answer.

% to utilize EK's code
%PsychCurve_general_ori(data(:,[6,5]), animalId)

%Sided data with outputs
 %[deltaAngle, pv, sl, pco, ci_thr, ci_s, ci_data] = PsychCurve_sidepref_ori(data(:,[6,5]) ,numSessions, animalId, 50)

%[StimLevel, NumPos, OutOfNum] = process_experiments_unsided(data(:,[6,5]))
 
% PsychCurve_sidepref_ori(experiments, numSessions,ferret,N)
% used in bootstrap to find confidence intervals

% %% Create the variables 'StimLevels', 'NumPos' and 'OutOfNum' for use in Palamedes Baysian fitting. 
% %This comes from the sided data to make a Logistic function. Called 'sided'
% %in the variables here
% 
% %unsided data across all trials
% cohLevels = round(data(:,4)*100); %making into intergers and rounding
% cohLevels=unique(cohLevels);
% correctTrial = data(:,5);
% 
% %Sided analysis
% cohLevelsSided=round((data(:,6)*100));
% cohLevelsSided=unique(cohLevelsSided);
% %StimLevels = cohLevelsSided';
% StimLevels = cohLevels';
% PosTs = data(:,7);  %did it go to the right
% 
% %NumPos is the number in which they responded by going to the positive spout (right sided answers  (-180 to 180)
% % Run this version for the variables of a sided curve
% % for i=1:length(StimLevels)
% %     idx=find((round(data(:,6)*100)) ==StimLevels(i));
% %     NumPos(i) = sum(PosTs(idx));
% % end
% % 
% % %OutOfNum is the total trials run at each coherence value
% % trialsRun = ones(length(correctTrial), 1);
% % for i=1:length(cohLevelsSided)
% %     idx=find((round(data(:,6)*100)) == cohLevelsSided(i));
% %    OutOfNum(i)=sum(trialsRun(idx)); 
% % end
% 
% 
% %For unsided analysis
% % Run this version for the variables of a Unsided (averaged) curve
% for i=1:length(cohLevels)
%     idx=find((round(data(:,4)*100)) ==cohLevels(i));
%     NumPos(i) = sum(data(idx,5));
% end
% 
% %OutOfNum is the total trials run at each coherence value
% trialsRun = ones(length(correctTrial), 1);
% for i=1:length(cohLevels)
%     idx=find((round(data(:,4)*100)) == cohLevels(i));
%    OutOfNum(i)=sum(trialsRun(idx)); 
% end
% 
% 
% %PF = @PAL_CumulativeNormal;
% %% summary
% %general evolution of coherence
% figure
% plot(data(:,4))
% xlabel('Coherence');
% ylabel('Trial Number');
% title(animalId, 'General Evolution of Coherence')
% 
% 
% %unsided data across all trials
% cohLevels = round(data(:,4)*100); %making into intergers and rounding
% cohLevels=unique(cohLevels); % for some reason "unique" gives duplicate rows for 0.05 to 0.20. try 
% correctTrial = data(:,5);
% 
% for i=1:length(cohLevels)
%     idx=find((round(data(:,4)*100)) == cohLevels(i));
%    propCorr(i)=sum(correctTrial(idx))/length(idx); 
%   
% end
% 
% 
% %Also want to add number of trials to plots to show some coherences are
% %rarely encounter and other are more dominant
% for i=1:length(cohLevels)
%     idx=find((round(data(:,4)*100)) == cohLevels(i));
%    numCorr(i)=sum(correctTrial(idx)); 
%   
% end
% 
% %orginal plot without scalling in accoradance with trial number
% % figure
% % 
% % plot( cohLevels, propCorr,'o') 
% % xlabel('Coherence');
% % ylabel('Percent of Trials');
% % title(animalId, 'Portion Correct')
% % %Note: E1 has gotten down to 0% coherence and guessed correctly several
% % %times which throws things off
% 
% 
% %Lets find the total number of trials at a given coherence
% 
% trialsRun = ones(length(correctTrial), 1);
% for i=1:length(cohLevels)
%     idx=find((round(data(:,4)*100)) == cohLevels(i));
%    trialNum(i)=sum(trialsRun(idx)); 
% end
% 
% figure
% %is a scatter plot with the total number of trials run at a particular
% %coherence scaling the size of the marker for each coherence
% scatter(cohLevels, (propCorr*100), trialNum)
% xlabel('Coherence');
% ylabel('Percent of Trials');
% title(animalId, 'Portion Correct')
% 
% 
% %Sided analysis
% cohLevelsSided=round((data(:,6)*100));
% cohLevelsSided=unique(cohLevelsSided);
% RsideTrial = data(:,7);
% 
% for i=1:length(cohLevelsSided)
%     idx=find((round(data(:,6)*100)) ==cohLevelsSided(i));
%     propCorrSide(i) = sum(RsideTrial(idx))/length(idx);
% end
% 
% for i=1:length(cohLevelsSided)
%     idx=find((round(data(:,6)*100)) == cohLevelsSided(i));
%    trialNumSided(i)=sum(trialsRun(idx)); 
% end
% %propCorrSide = propCorrSide.'; %i dont know why it comes out going the wrong way but this fixes it
% 
% 
% figure
% %plots to both negative and positive sides
%  scatter(cohLevelsSided, (propCorrSide*100), trialNumSided)
% xlabel('Coherence');
% ylabel('Percent of Trials');
% title(animalId, 'Trials to Side 1 by Direction')
% 
% 
% % plot( cohLevelsSided, propCorrSide,'o')
% % xlabel('Coherence');
% % ylabel('Portion of Trials to Spout 1');
% % title(animalId, 'Portion Correct by Side')
% 
% 
% %% Save data
% %
% sessionResults.files = fnames;
% sessionResults.data = data;
% sessionResults.cohLevels = cohLevels;
% sessionResults.cohLevelsSided = cohLevelsSided;
% sessionResults.propCorrect = propCorr;
% sessionResults.propCorrectSided = propCorrSide;
% sessionResults.PDS = PDS;
% 
% cd 'Z:\Ferret Behavior\RDK\Amblyopia\stairs dominant'
% save(animalId,'sessionResults')
% 
