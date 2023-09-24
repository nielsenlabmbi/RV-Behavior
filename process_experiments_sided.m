function [StimLevels, NumRight, OutOfNum] = process_experiments_sided(experiments)

stim = experiments(:,1);
stimLeft = stim(stim < 0);
stimRight = stim(stim > 0);

StimLevelsLeft =  unique(stimLeft)';
StimLevelsRight =  unique(stimRight)';
StimLevels = [StimLevelsLeft StimLevelsRight];

%Number of positive responses (e.g., 'yes' or 'correct' at each of the 
%   entries of 'StimLevels'  
%NumPos = [60 77.7778 50 85.7143 87.5 88.8889 87.5 88.8889 81.8182 100 100 100];                 
NumPosLeft = [];
for i=1:length(StimLevelsLeft);
    NumPosLeft = [NumPosLeft length(find(experiments(ismember(stim,StimLevelsLeft(i)),2)))];
end

NumPosRight = [];
for i=1:length(StimLevelsLeft);
    NumPosRight = [NumPosRight length(find(experiments(ismember(stim,StimLevelsRight(i)),2)))];
end
%Number of trials at each entry of 'StimLevels'
%OutOfNum = [10 9 2 7 8 9 8 9 11 6 10 11];

OutOfNumLeft = [];
for i=1:length(StimLevelsLeft);
    OutOfNumLeft = [OutOfNumLeft length(find(ismember(stim,StimLevelsLeft(i))))];
end

NumNegLeft = OutOfNumLeft - NumPosLeft;

NumRight = [NumNegLeft NumPosRight];

OutOfNum = [];
for i=1:length(StimLevels);
    OutOfNum = [OutOfNum length(find(ismember(stim,StimLevels(i))))];
end

end