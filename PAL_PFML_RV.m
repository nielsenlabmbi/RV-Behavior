%
%PAL_PFML_Demo  Demonstrates use of Palamedes functions to (1) fit a
%Psychometric Function to some data using a Maximum Likelihood criterion, 
%(2) determine standard errors of free parameters using a bootstrap
%procedure and (3) determine the goodness-of-fit of the fit.
%
%Demonstrates basic usage of Palamedes functions:
%-PAL_PFML_Fit
%-PAL_PFML_BootstrapParametric
%-PAL_PFML_BootstrapNonParametric
%-PAL_PFML_GoodnessOfFit
%secondary:
%-PAL_Logistic
%
%More information on any of these functions may be found by typing
%help followed by the name of the function. e.g., help PAL_PFML_Fit
%
%NP (08/23/2009)

tic

message = 'Parametric Bootstrap (1) or Non-Parametric Bootstrap? (2): ';
ParOrNonPar = input(message);

% Example input 
% %Stimulus intensities
% StimLevels = [0.01 0.03 0.05 0.07 0.09 0.11]; 
% 
% %Number of positive responses (e.g., 'yes' or 'correct' at each of the 
% %   entries of 'StimLevels'  
% NumPos = [59 53 68 83 92 99];                 
% 
% %Number of trials at each entry of 'StimLevels'
% OutOfNum = [100 100 100 100 100 100];         

%Use the Logistic function
PF =@PAL_CumulativeNormal ;  %Alternatives: @PAL_Gumbel, @PAL_Weibull,
                     %@PAL_Quick, @PAL_logQuick, @PAL_Logistic
                     %@PAL_CumulativeNormal, @PAL_HyperbolicSecant

%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
 
%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
%*Here the original demo parameters
% searchGrid.alpha = 0.01:.001:.11;
% searchGrid.beta = logspace(0,3,101);
% searchGrid.gamma = 0.15;  %scalar here (since fixed) but may be vector
% searchGrid.lambda = 0.02;  %ditto

%****Here using the values for Erika's dots figure; CHANGE LATER***
searchGrid.alpha = 0.2:0.002:1; %threshold Is this what somehow gets set to 85 and 72? Like it would end at .85?
searchGrid.beta = logspace(-1,1,200); %slope
searchGrid.gamma = [0, 0.1, 1];  %RV made vector because curve just looks wrong and has 0.5 as an asymtope when scalar 
% and it is 0.5 or 0.2 everywhere else I see it.
searchGrid.lambda = floor(100-max(NumPos(end)/OutOfNum(end))*100)/100;  %can also be vector
%%
%Perform fit 
%LL things are 'Log likelyhoods'
disp('Fitting function.....');
[paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
    OutOfNum,searchGrid,paramsFree,PF);

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValues(1));  % the %6.4f means when prints out the result it will have 6 digits before the decimal, a decimal, and 4 digits after.
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
disp(message);

%Number of simulations to perform to determine standard error
B=400;                  

disp('Determining standard errors.....');

if ParOrNonPar == 1
    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(... 
        StimLevels, OutOfNum, paramsValues, paramsFree, B, PF, ...
        'searchGrid', searchGrid);
else
    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
        StimLevels, NumPos, OutOfNum, [], paramsFree, B, PF,...
        'searchGrid',searchGrid);
end

disp('done:');
message = sprintf('Standard error of Threshold: %6.4f',SD(1));
disp(message);
message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
disp(message);

%Distribution of estimated slope parameters for simulations will be skewed
%(type: hist(paramsSim(:,2),40) to see this). However, distribution of
%log-transformed slope estimates will be approximately symmetric
%[type: hist(log10(paramsSim(:,2),40)]. This might motivate using 
%log-scale for slope values (uncomment next three lines to put on screen):
% SElog10slope = std(log10(paramsSim(:,2)));
% message = ['Estimate for log10(slope): ' num2str(log10(paramsValues(2))) ' +/- ' num2str(SElog10slope)];
% disp(message);
%%
%Number of simulations to perform to determine Goodness-of-Fit
B=1000;

disp('Determining Goodness-of-fit.....');

[Dev pDev] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
    paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);

disp('done:');

%Put summary of results on screen
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('p-value: %6.4f',pDev);
disp(message);
 
%Create simple plot
ProportionCorrectObserved=NumPos./OutOfNum; 
StimLevelsFineGrain=[min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
 
figure('name','Maximum Likelihood Psychometric Function Fitting');
axes
hold on
plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[0 .7 0],'linewidth',4);
plot(StimLevels,ProportionCorrectObserved,'k.','markersize',40);
set(gca, 'fontsize',16);
set(gca, 'Xtick',StimLevels);
axis([min(StimLevels) max(StimLevels) 0 1]);
xlabel('Dot Coherence');
ylabel('Fraction to Right');


%Now trying to figure out what is meant by having threshold evaluations
%  " x = PAL_Weibull(params, y, 'Inverse') returns the x-value at 
%   which the Psychometric Function evaluates to y." So below finds the
%   threshold giving .75 meaning 75% accuarate, may also be 75 depending on
%   how the data is entered.

Threshold  = PAL_Weibull(paramsValues,.75,'Inverse');
toc