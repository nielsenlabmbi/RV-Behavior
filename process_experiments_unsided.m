function [StimLevels, NumPos, OutOfNum] = process_experiments_unsided(experiments)

StimLevels =  unique(abs(experiments(:,1)))';
% Number of positive responses (e.g., 'yes' or 'correct' at each of the 
%   entries of 'StimLevels'  
       
NumPos = [];
for i=1:length(StimLevels);
    NumPos = [NumPos length(find(experiments(ismember(abs(experiments(:,1)),StimLevels(i)),2)))];
end
% Number of trials at each entry of 'StimLevels'

OutOfNum = [];
for i=1:length(StimLevels);
    OutOfNum = [OutOfNum length(find(ismember(abs(experiments(:,1)),StimLevels(i))))];
end

end