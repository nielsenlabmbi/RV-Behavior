%function [Threshold, StimLevelsFineGrain,ProportionCorrectModel, StimLevels, ProportionCorrectObserved, PF, paramsValues] = PsychCurve_general_ori(OutOfNum, NumPos, StimLevels, animalId)
function [Threshold, StimLevelsFineGrain,ProportionCorrectModel, StimLevels, ProportionCorrectObserved, PF, paramsValues] = PsychCurve_sided_RDK(experiments, animalId)
%%%%% PALAMEDES FUNCTIONS FOR BOTH TOGETHER%%%%%%%%
%is making figure 3.d from the EK 2019 paper.   This is not a sigmoid curve
%but the saturation curve looking kagigger.
%this is grabbing the averaged across left or right or correct to just show
%the correct at a given stimulus

%EK original experiments matrix is a structure of results which can be used
%to calculate the parameters. In versions as of 9/10/23 these are generated
%in a seperate script then called 

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
%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = 0:0.1:30; %threshold
searchGrid.beta = logspace(0,0.1,100); %slope
searchGrid.gamma = .50;  %guess rate: is set to 0.5 because Left vs Right is a binary choice. scalar here (since fixed) but may be vector
searchGrid.lambda = 0.02;  %lapse rate scalar when fixed fixed, vector when free

%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
 
%Fit a Logistic function
PF = @PAL_Weibull;  %Alternatives: PAL_Gumbel, PAL_Weibull, 
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant

%Optional:  RV What does this do?
% these variables(?terminology?) are present in the lapseFit_Demo.  Is that
% what is happening here?
options = PAL_minimize('options');   %type PAL_minimize('options','help') for help
options.TolFun = 1e-09;     %increase required precision on LL = tolerance for changes in function.
options.MaxIter = 100;     % max iteration number, limits number of times it can run
options.Display = 'off';    %suppress fminsearch messages
%fminsearch finds arguments that give the minimum value for parameters of a
%function

%Perform fit
disp('Fitting function.....');
[paramsValues, LL, exitflag, output] = PAL_PFML_Fit(StimLevels,NumPos, ...
    OutOfNum,searchGrid,paramsFree,PF,'searchOptions',options);

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
disp(message);




%Create simple plot %does not include confidence intervals
ProportionCorrectObserved=NumPos./OutOfNum; 
StimLevelsFineGrain=(0:max(StimLevels)/1000:max(StimLevels));
ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);

%get the threshold to reach 75% correct, change to .82 for 82% correct
Threshold  = PAL_Weibull(paramsValues,.75,'Inverse');

figure('name','Maximum Likelihood Psychometric Function Fitting');
plot(StimLevels,ProportionCorrectObserved,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',StimLevels);
axis([0 max(StimLevels) 0 1]);
hold on;
plot(StimLevelsFineGrain,ProportionCorrectModel,'g-','linewidth',4);
hold on; plot(StimLevels,.5*(ones(size(StimLevels))), 'r--') %just makes a red line at 50%
xlabel('coherence level');
ylabel('proportion correct');
set(gca, 'fontsize',11);
legend (num2str(Threshold),'Location','southeast')
title(strcat(animalId,'  ', num2str(round(mean(OutOfNum))),...
    ' trials per stimulus'))


