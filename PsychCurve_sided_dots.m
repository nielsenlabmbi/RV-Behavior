
function [paramsValues,deltaCoherence, Dev, pDev, ci_thr, ci_s, ci_data, data_boot, t_boot, sl_boot] =  PsychCurve_sided_dots(StimLevels, NumPos, OutOfNum,animalId,N,pDevBool)
%%%%% PALAMEDES FUNCTIONS FOR %to the Right /SIDED ANALYSIS%%%%%%%%

%Required inputs: 
% StimLevels: absolute value of conherences 
% NumPos: number correct at each StimLevel 
% OutOfNum: total number of trials run
% animalId,
% N: number of iterations for the bootstrap to obtain confidence intervals 
% pDevBool: include as a variable to run Goodness of Fit

% Generally assume to do Non-parametric as biases will vary and may prevent
% Gaussian distribution.
message = 'Parametric Bootstrap (1) or Non-Parametric Bootstrap? (2): ';
ParOrNonPar = input(message);

%Stimulus intensities

% From EK's version where data is set in a structure called 'experiments'
% rather than created as arguments durring initial processing.

% StimLevels =  unique(experiments(:,1))';

% %Number of positive responses (e.g., 'yes' or 'correct' at each of the 
% %   entries of 'StimLevels'               
% NumPos = [];
% for i=1:length(StimLevels);
%     NumPos = [NumPos length(find(experiments(ismember(experiments(:,1),StimLevels(i)),2)))];
% end

% %Number of trials at each entry of 'StimLevels'
% OutOfNum = [];
% for i=1:length(StimLevels);
%     OutOfNum = [OutOfNum length(find(ismember(experiments(:,1),StimLevels(i))))];
% end

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = 0:0.002:1;  %search through these for alpha, the threshold
searchGrid.beta = logspace(-1,1,200); % search through here for slope
% gamma can be lower if ferret has bias that makes performance worse than
% chance (gamma = 0.50)
%searchGrid.gamma = 0.3:0.05:0.6; %for when possibly doing worse than just guessing
searchGrid.gamma = 0.0:0.05:0.3; %this sets left lapse rate to scan lowest possible range
%searchGrid.lambda = floor(100-max(NumPos(end)/OutOfNum(end))*100)/100; %this sets right side lapse rate to highest % Correct
%searchGrid.lambda = 0.0:0.05:0.3; %this sets left lapse rate to scan lowest possible range
searchGrid.lambda = 0.3:0.10:0.7; %range for when worse than guess, eg. lesion studies
%Alternative parameters
% searchGrid.alpha = 1:.1:30;
%searchGrid.beta = logspace(-1,3,500);
%searchGrid.gamma = 0.5; % use for unsided as 50% represents rate if only guessing.  
% searchGrid.gamma = NumPos(1)/OutOfNum(1);%0.25:.05:.75;%.50;  %scalar here (since fixed) but may be vector
% searchGrid.lambda = 1 - NumPos(end)/OutOfNum(end); %%0:.05:.5;%0.02;  %ditto

%Threshold and Slope and lapse rates are fixed
paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter
 
%Fit a Logistic function
PF = @PAL_CumulativeNormal;     %Alternatives:PAL_Weibull, PAL_Gumbel, PAL_Logistic, 
                     %PAL_CumulativeNormal, PAL_HyperbolicSecant, but EK
                     %used Weibull so we will stick with that

%Optional:
options = PAL_minimize('options');   %type PAL_minimize('options','help') for help
options.TolFun = 1e-09;     %increase required precision on LL
options.MaxIter = 100;
options.Display = 'off';    %suppress fminsearch messages

%Perform fit
disp('Fitting function.....');
[paramsValues, paramsValues(5), paramsValues(6)]  = PAL_PFML_Fit(StimLevels,NumPos, ...
    OutOfNum,searchGrid,paramsFree,PF,'searchOptions',options,'LapseLimits',[0 1],'GuessLimits',[0 1]);

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
disp(message);

%Bootstrap confidence intervals
%if exist('N','var') & N > 0
[ci_thr, ci_s, ci_data, data_boot, t_boot, sl_boot] = bootstrap_pf(StimLevels,NumPos, OutOfNum, N, searchGrid, paramsFree,PF,options);
%end

% Used in unsided % correct calculations
% %Find Threhold for 75% correct
% % Threshold75  = PAL_Weibull(paramsValues,.75,'Inverse');
% % message = ['Threshold estimate at 75%: ' num2str(Threshold75)];
% % disp(message);
% 
% %Find Threhold for 82% correct
% % Threshold82  = PAL_Weibull(paramsValues,.82,'Inverse');
% % message = ['Threshold estimate at 82%: ' num2str(Threshold82)];
% % disp(message);

% Find coherence needed to shift from 50-68% correct.
upperbound = 1-paramsValues(4); %paramsValues(4) is lapse rate, thus 1-paramsValues(4) is the upperbound on performance. 
guessrate = paramsValues(3); %percentage at 0
%for sided data use CumulativeNormal, for unsided use Weibull

middlepoint = PAL_CumulativeNormal(paramsValues,guessrate+(upperbound-guessrate)/2,'Inverse'); 
sixtyseven = PAL_CumulativeNormal(paramsValues,guessrate+2*(upperbound-guessrate)/3,'Inverse');
deltaCoherence = sixtyseven - middlepoint;

% %Number of simulations to perform to determine standard error
% B=400;                  
% 
% disp('Determining standard errors.....');
% 
% if ParOrNonPar == 1
%     [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(...
%         StimLevels, OutOfNum, paramsValues, paramsFree, B, PF, ...
%         'searchOptions',options,'searchGrid', searchGrid);
% else
%     [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
%         StimLevels, NumPos, OutOfNum, [], paramsFree, B, PF,...
%         'searchOptions',options,'searchGrid',searchGrid);
% end
% 
% disp('done:');
% message = sprintf('Standard error of Threshold: %6.4f',SD(1));
% disp(message);
% message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
% disp(message);
% 
%
if exist('pDevBool','var') & pDevBool
%Number of simulations to perform to determine Goodness-of-Fit
B=400;
% 
disp('Determining Goodness-of-fit.....');
%Output: 
%   'Dev': Deviance (transformed likelihood ratio comparing fit of
%       psychometric function to fit of saturated model)
%
%   'pDev': proportion of the B Deviance values from simulations that were
%       greater than Deviance value of data. The greater the value of pDev,
%       the better the fit.

[Dev, pDev] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
    paramsValues(1:4), paramsFree, B, PF,'searchOptions',options, ...
    'searchGrid', searchGrid)
% 
disp('done:');

%Put summary of results on screen
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('p-value: %6.4f',pDev);
disp(message);

else
    pDev = NaN; 
end
% %Create simple plot
ProportionCorrectObserved=NumPos./OutOfNum; 
StimLevelsFineGrain=(min(StimLevels):max(StimLevels)/1000:max(StimLevels));
ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
%  
% figure('name','Maximum Likelihood Psychometric Function Fitting');
% plot(StimLevels,ProportionCorrectObserved,'k.','markersize',40);
% set(gca, 'fontsize',16);
% set(gca, 'Xtick',StimLevels);
% axis([min(StimLevels) max(StimLevels) 0 1]);
% hold on;
% plot(StimLevelsFineGrain,ProportionCorrectModel,'g-','linewidth',4);
% hold on; plot(StimLevels,.5*(ones(size(StimLevels))), 'r--')
% errorbar(StimLevels, ProportionCorrectObserved, ((ci_data(2,:)) - (ci_data(1,:))), ((ci_data(3, :)) - (ci_data(2,:))))
% xlabel('coherence level');
% ylabel('proportion correct');
% set(gca, 'fontsize',11);
% %boot_data is the 95% confidence intervals
%title(animalId)

