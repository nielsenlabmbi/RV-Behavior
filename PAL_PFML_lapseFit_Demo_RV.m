%
%PAL_PFML_lapseFit_Demo demonstrates the options for incorporating the 
%lapse rate into a PF model. These options were introduced in Palamedes 
%version 1.3.0.

%RV Needs to figure out what portions are demo and what are part of the
%fitting process.   Then what needs to be called in the lab's code to make the
%fitting happen
%
%Demonstrates usage of Palamedes functions:
%-PAL_PFML_Fit
%-PAL_PFML_BootstrapParametric
%-PAL_PFML_BootstrapNonParametric
%-PAL_PFML_GoodnessOfFit
%-PAL_PFML_FitMultiple
%-PAL_PFML_BootstrapParametricMultiple
%-PAL_PFML_BootstrapNonParametricMultiple
%-PAL_PFML_GoodnessOfFitMultiple
%-PAL_PFLR_ModelComparison
%-PAL_PF_SimulateObserverParametric
%The above are demonstrated here especially with regard to usage of the
%optional function argument 'lapseFit'.
%
%NP (September 2011)

%clear all

rand('state',sum(100*clock));

%Nelder-Mead search options
options = PAL_minimize('options');  %decrease tolerance (i.e., increase
options.TolX = 1e-09;              %precision).
options.TolFun = 1e-09;  %RV no idea what these variables are. Does 'Fun' stand for function? 
options.MaxIter = 10000;  
options.MaxFunEvals = 10000;
%
question = sprintf('\nSingle or Multiple condition [S/M]?: ');
wish = input(question,'s');
% RV I think this means single or multiple variables in the experiment. We have side and coherence so M is the answer? 

if strncmpi(wish,'S',1)

    disp([sprintf('\n') 'A psychometric Function will be fitted to some simulated data.' sprintf('\n')]);
    disp(['Lapse rates can be handled in any of 3 ways. ''nAPLE'' jointly fits thresholds,']);
    disp(['slopes, and lapse rates, in the manner advocated by Wichmann & Hill (2001a).']);
    disp(['''jAPLE'' is identical to ''nAPLE'' except that it assumes that the highest']);
    disp(['stimulus intensity is at an asymptotic level and thus that an error observed']);
    disp(['at this intensity can only be due to lapses.' sprintf('\n')]);    
 %RV not stated in lab documents which is done but jAPLE  sounds most relevent to our observations so we can do that.   
    lapseFit = '';
    while ~strncmpi(lapseFit,'nAPLE',3) && ~strncmpi(lapseFit,'jAPLE',3)
        lapseFit = input('How do you wish to fit lapse rates? [nAPLE/jAPLE]?: ','s');
    end
%
    %prepare plot %RV but what is it plotting? No data has been entered?
    figure('units','pixels','position',[100 100 400 200]);
    h = axes('units','pixels','position',[50 50 300 100]);
    set(h,'xtick',[-2:1:2 10],'xlim',[-3 11],'ylim',[.4 1],'xticklabel',{'-2','-1','0','1','2','10'});
    xlabel('Stimulus Intensity');    
    ylabel('prop correct');    
    hold on;  %Ah its just preping to make a figure with data it will run. This is backward from how I would do it 
    %and the parameters should exist before the fitting
    
    
    
    %Stimulus intensities. Generating logistic (F) evaluates to near unity
    %at 10
   % makes data for running the demo.
    StimLevels = [-2 -1 0 1 2 10];

    OutOfNum = [150 150 150 150 150 150];    %N need not be equal

    PF = @PAL_CumulativeNormal;                     %PF function
    paramsValues = [0 1 .5 .05];            %generating values

    %Simulate observer  %RV - the words in '-' just say to look for the
    %parammeter afterthem not in '-'.
    
    %NumPos is also something I load.
    NumPos = PAL_PF_SimulateObserverParametric(paramsValues,StimLevels,OutOfNum,PF,'lapseFit',lapseFit);
    plot(h,StimLevels,NumPos./OutOfNum,'ko','markersize',6,'markerfacecolor','k');
    
    %Fit PF
    %More toy information i think. It does not need these particular values
    %to run.
    searchGrid.alpha = [-1:.05:1];    %structure defining grid to
    searchGrid.beta = 10.^[-1:.05:2]; %search for initial values
    searchGrid.gamma = .5;
    searchGrid.lambda = [0:.005:.1];

    paramsFree = [1 1 0 1]; %[threshold slope guess lapse] 1: free, 0:fixed
    
    [paramsFitted LL exitflag] = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, searchGrid, paramsFree, PF,'lapseLimits',[0 1],'searchOptions',options,'lapseFit',lapseFit);

    if ~exitflag
        disp('Psychometric Function Fit failed! Exiting...');
        return
    end
    
    plot(h,-2.5:.01:10.5,PF(paramsFitted,-2.5:.01:10.5),'-','color',[0 .7 0],'linewidth',2);
    
    disp(sprintf('\n'));   
    disp(['Threshold: ' num2str(paramsFitted(1),'%4.3f')]);
    disp(['Slope: ' num2str(paramsFitted(2),'%4.3f')]);
    disp(['Guess Rate: ' num2str(paramsFitted(3),'%4.3f')]);
    disp(['Lapse Rate: ' num2str(paramsFitted(4),'%4.3f')]);
    disp(sprintf('\n'));    

    
    question = sprintf('\nDo you wish to perform a non-parametric bootstrap to obtain\nSEs for parameters [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
    
        [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(StimLevels, NumPos, OutOfNum, ...
        [],paramsFree, B, PF, 'lapseLimits',[0 1], 'searchGrid', searchGrid,'lapseFit',lapseFit);
    
        disp([ sprintf('\n') 'SE for the threshold: ' num2str(SD(1),'%4.3f')]);
        disp(['SE for the slope: ' num2str(SD(2),'%4.3f')]);
        disp(['SE for the lapse rate: ' num2str(SD(4),'%4.4f') sprintf('\n')]);
    end
    
    question = sprintf('\nDo you wish to perform a parametric bootstrap to obtain\nSEs for parameters [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
    
        [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(StimLevels, OutOfNum, ...
        paramsFitted, paramsFree, B, PF, 'lapseLimits',[0 1], 'searchGrid', searchGrid,'lapseFit',lapseFit);
    
        disp([ sprintf('\n') 'SE for the threshold: ' num2str(SD(1),'%4.3f')]);
        disp(['SE for the slope: ' num2str(SD(2),'%4.3f')]);
        disp(['SE for the lapse rate: ' num2str(SD(4),'%4.4f') sprintf('\n')]);
    end

    %would this be fore the whole curve or just the lapse rates? It looks
    %like the whole curve
    question = sprintf('\nDo you wish to perform a Goodness-of-fit test [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
    
        [Dev pDev DevSim converged] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
        paramsFitted, paramsFree, B, PF, 'searchGrid', searchGrid, 'lapseLimits',[0 1],'lapseFit',lapseFit);
      
        disp([sprintf('\n') 'Goodness-of-fit by Monte Carlo: ' num2str(pDev,'%4.4f')])

        if exist('chi2pdf.m') == 2
            disp(['Goodness-of-fit by chi-square approximation: ' num2str(1-chi2cdf(Dev,6-sum(paramsFree)),'%4.4f')])
        end
    end
    
end

%Can I trials to the left and trials to the right as two differnt
%conditions like
% leftStim = 
%if data(:,4) >=0
if strncmpi(wish,'M',1)

    disp([sprintf('\n') 'Two models will be fit to PFs obtained in 4 conditions. The fuller']);
    disp(['model allows each condition to have its own threshold and slope, but assumes']);
    disp(['lapse rates are equal. The lesser model assumes that threshold, slopes, and']);
    disp(['lapse rates are equal across conditions.' sprintf('\n')]);    
    disp(['Lapse rates can be handled in any of 3 ways. ''nAPLE'' jointly fits thresholds,']);
    disp(['slopes, and lapse rates, in the manner advocated by Wichmann & Hill (2001a).']);
    disp(['''jAPLE'' is identical to ''nAPLE'' except that it assumes that the highest']);
    disp(['stimulus intensity is at an asymptotic level and thus that an error observed']);
    disp(['at this intensity can only be due to lapses.' sprintf('\n')]);    
    % none of these sound like Erika's "seperate lapse rates were fit for
    % left and right responses". reading I guess we do jAPLE? ferrets could
    % get 100% at 100% coherence.

    lapseFit = '';
    while ~strncmpi(lapseFit,'nAPLE',3) && ~strncmpi(lapseFit,'jAPLE',3)
        lapseFit = input('How do you wish to fit lapse rates? [nAPLE/jAPLE]?: ','s');
    end

    %prepare plots   --RV do not see where "cond" is made.
    figure('units','pixels','position',[100 100 400 450]);
    for cond = 1:4
        h(cond) = axes('units','pixels','position',[50 350-(cond-1)*100 300 75]);%125 300 100]);
        set(h(cond),'xtick',[-2:1:2 10],'xlim',[-3 11],'ylim',[.4 1],'xticklabel',{'-2','-1','0','1','2','10'});
        ylabel('prop correct');    
        hold on;
    end
    
    xlabel('Stimulus Intensity');    
    h(5) = axes('units','pixels','position',[50 425 300 40]);
    set(h(5),'xlim',[0 1],'ylim',[0 1]);
    line([.05 .15],[.5 .5],'color',[.7 0 0],'linewidth',2)
    text(.2,.5,'lesser model');
    line([.5 .6],[.5 .5],'color',[0 .7 0],'linewidth',2)
    text(.65,.5,'fuller model');
    axis off;

    %Stimulus intensities. Generating logistic (F) evaluates to near unity
    %at 10
 % toy data for running as a demo i think   
    StimLevels = repmat([-2 -1 0 1 2 10],[4 1]);

    OutOfNum = [100 100 100 100 100 100;    %N need not be equal
               150 150 150 150 150 150;
               100 100 100 100 100 100;
               150 150 150 150 150 150];

    PF = @PAL_Logistic;                     %PF function
    paramsValues = [0 1 .5 .05];            %generating values

    %Simulate observer  Again- is this par of the analysis or the demo? Do
    %i need to simulate and observer for compairison or will that wipe out the actual
    %data?
    
    %Code is written to model multiple conditions?  Are these like multiple
    %experiments? Eg. Cond1 is dotSpeed 72 and cond2 is dotSpeed 48? In the
    %model all of the condtiosns are the same 
    % OO or based on looking the the NumPos are these broken down by an
    % equivelant of Sessions? So I would feed in NumPos and all as
    % matriceis? No I think it can be run on just arrays if looking at just
    % one animal over multiple sessions in one condition. NumPos and all
    % would then just be an array and some stuff may need to be tweaked.
    for cond = 1:4
        NumPos(cond,:) = PAL_PF_SimulateObserverParametric(paramsValues,StimLevels(cond,:),OutOfNum(cond,:),PF,...
            'lapseFit',lapseFit);
        plot(h(cond),StimLevels(cond,:),NumPos(cond,:)./OutOfNum(cond,:),'ko','markersize',6,'markerfacecolor','k');
    end

    %Define fuller model
    thresholdsfuller = 'unconstrained';  %Each condition gets own threshold
    slopesfuller = 'unconstrained';      %Each condition gets own slope
    guessratesfuller = 'fixed';          %Guess rate fixed, RV is there a way to make it not fixed? What are the phrases for that?
    lapseratesfuller = 'unconstrained';    %'constrained' means the lapse rates are the same, 
    %'unconstrained' calculates the sides of the cureve independently

    %Fit fuller model
    [paramsFuller LL exitflag trash trash numParamsFuller] = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
      paramsValues, PF,'searchOptions',options,'lapserates',lapseratesfuller,'thresholds',thresholdsfuller,...
      'slopes',slopesfuller,'guessrates',guessratesfuller,'lapseLimits',[0 1],'lapseFit',lapseFit);
  
    if ~exitflag
        disp('Fit to fuller model failed! Exiting...');
        return
    end
    
    
    disp(sprintf('\n'))
    disp('Fuller Model:')
    disp(sprintf('\n'));
    disp(['Thresholds: ' num2str(paramsFuller(1,1),'%4.3f') ', ' num2str(paramsFuller(2,1),'%4.3f') ', ' num2str(paramsFuller(3,1),'%4.3f') ', ' num2str(paramsFuller(4,1),'%4.3f')]);
    disp(['Slopes: ' num2str(paramsFuller(1,2),'%4.3f') ', ' num2str(paramsFuller(2,2),'%4.3f') ', ' num2str(paramsFuller(3,2),'%4.3f') ', ' num2str(paramsFuller(4,2),'%4.3f')]);
    disp(['Guess Rates: ' num2str(paramsFuller(1,3),'%4.3f') ', ' num2str(paramsFuller(2,3),'%4.3f') ', ' num2str(paramsFuller(3,3),'%4.3f') ', ' num2str(paramsFuller(4,3),'%4.3f')]);
    disp(['Lapse Rates: ' num2str(paramsFuller(1,4),'%4.3f') ', ' num2str(paramsFuller(2,4),'%4.3f') ', ' num2str(paramsFuller(3,4),'%4.3f') ', ' num2str(paramsFuller(4,4),'%4.3f')]);
    disp(sprintf('\n'));
    disp(['Akaike''s Informaton Criterion: ' num2str(-2*LL + 2*numParamsFuller,'%4.3f')]);
    disp(sprintf('\n'));
    %it has all the params in plural, why?  Akaike's Information Cri is a
    %way select a model. The higher the value, the less information is lost
    %in the model
    
    %plot fuller model
    for cond = 1:4
        plot(h(cond),-2.5:.01:10.5,PF(paramsFuller(cond,:),-2.5:.01:10.5),'-','color',[0 .7 0],'linewidth',2);
    end
    drawnow

    %Define lesser model
    thresholdslesser = 'constrained';   %Common threshold
    slopeslesser = 'constrained';       %Common slope
    guessrateslesser = 'fixed';         %Guess rate fixed
    lapserateslesser = 'unconstrained';   %Common lapse rate

    %Fit lesser model
    [paramsLesser LL exitflag trash trash numParamsLesser] = PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, ...
      paramsValues, PF,'searchOptions',options,'lapserates',lapserateslesser,'thresholds',thresholdslesser,...
      'slopes',slopeslesser,'guessrates',guessrateslesser,'lapseLimits',[0 1],'lapseFit',lapseFit) ;

    if ~exitflag
        disp('Fit to lesser model failed! Exiting...');
        return
    end
    
    disp(sprintf('\n'))
    disp('Lesser Model:')
    disp(sprintf('\n'));
    disp(['Thresholds: ' num2str(paramsLesser(1,1),'%4.3f') ', ' num2str(paramsLesser(2,1),'%4.3f') ', ' num2str(paramsLesser(3,1),'%4.3f') ', ' num2str(paramsLesser(4,1),'%4.3f')]);
    disp(['Slopes: ' num2str(paramsLesser(1,2),'%4.3f') ', ' num2str(paramsLesser(2,2),'%4.3f') ', ' num2str(paramsLesser(3,2),'%4.3f') ', ' num2str(paramsLesser(4,2),'%4.3f')]);
    disp(['Guess Rates: ' num2str(paramsLesser(1,3),'%4.3f') ', ' num2str(paramsLesser(2,3),'%4.3f') ', ' num2str(paramsLesser(3,3),'%4.3f') ', ' num2str(paramsLesser(4,3),'%4.3f')]);
    disp(['Lapse Rates: ' num2str(paramsLesser(1,4),'%4.3f') ', ' num2str(paramsLesser(2,4),'%4.3f') ', ' num2str(paramsLesser(3,4),'%4.3f') ', ' num2str(paramsLesser(4,4),'%4.3f')  sprintf('\n')]);
    disp(['Akaike''s Informaton Criterion: ' num2str(-2*LL + 2*numParamsLesser,'%4.3f')  sprintf('\n')]);
    %okay, here is will display "Lapse Rates" that is data I want to
    %capture so I need to make an array or matrix of the paramsLesser and
    %paramsFuller(:,4)
    
    
    %plot lesser model
    for cond = 1:4
        plot(h(cond),-2.5:.01:10.5,PF(paramsLesser(cond,:),-2.5:.01:10.5),'-','color',[.7 0 0],'linewidth',2);
    end
    drawnow

    question = sprintf('\nDo you wish to perform a non-parametric bootstrap to obtain\nSEs for parameters of, say, the fuller model [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
        [SD paramsSim LLSim converged SDfunc funcParamsSim] = PAL_PFML_BootstrapNonParametricMultiple(StimLevels, ...
            NumPos, OutOfNum, paramsFuller, B, PF,'searchOptions',options,'lapserates',lapseratesfuller,'thresholds',...
            thresholdsfuller,'slopes',slopesfuller,'guessrates',guessratesfuller,'lapseLimits',[0 1],'lapseFit',lapseFit);
    
        disp([ sprintf('\n') 'SEs for the 4 thresholds: ' num2str(SD(1,1),'%4.3f') ', ' num2str(SD(2,1),'%4.3f') ', ' num2str(SD(3,1),'%4.3f') ', ' num2str(SD(4,1),'%4.3f')])
        disp(['SEs for the 4 slopes: ' num2str(SD(1,2),'%4.3f') ', ' num2str(SD(2,2),'%4.3f') ', ' num2str(SD(3,2),'%4.3f') ', ' num2str(SD(4,2),'%4.3f')])
        disp(['SE for the shared lapse rate: ' num2str(SD(1,4),'%4.4f')])
    end
    
    question = sprintf('\nDo you wish to perform a parametric bootstrap to obtain\nSEs for parameters of, say, the fuller model [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
        [SD paramsSim LLSim converged SDfunc funcParamsSim] = PAL_PFML_BootstrapParametricMultiple(StimLevels, ...
            OutOfNum, paramsFuller, B, PF,'searchOptions',options,'lapserates',lapseratesfuller,'thresholds',...
            thresholdsfuller,'slopes',slopesfuller,'guessrates',guessratesfuller,'lapseLimits',[0 1],'lapseFit',lapseFit);
        
        disp([sprintf('\n') 'SEs for the 4 thresholds: ' num2str(SD(1,1),'%4.3f') ', ' num2str(SD(2,1),'%4.3f') ', ' num2str(SD(3,1),'%4.3f') ', ' num2str(SD(4,1),'%4.3f')])
        disp(['SEs for the 4 slopes: ' num2str(SD(1,2),'%4.3f') ', ' num2str(SD(2,2),'%4.3f') ', ' num2str(SD(3,2),'%4.3f') ', ' num2str(SD(4,2),'%4.3f')])
        disp(['SE for the shared lapse rate: ' num2str(SD(1,4),'%4.4f')])
    end

    question = sprintf('\nDo you wish to determine Goodness-of-fit for, say, the fuller model [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
        [Dev, pDev DevSim converged] = PAL_PFML_GoodnessOfFitMultiple(StimLevels, NumPos, OutOfNum,paramsFuller, ...
            B, PF, 'searchOptions',options,'lapserates',lapseratesfuller,'thresholds',thresholdsfuller,'slopes',...
            slopesfuller,'guessrates',guessratesfuller,'lapseLimits',[0 1],'lapseFit',lapseFit);
        
        disp([sprintf('\n') 'Goodness-of-fit by Monte Carlo: ' num2str(pDev,'%4.4f')])
        if exist('chi2pdf.m') == 2
            disp(['Goodness-of-fit by chi-square approximation: ' num2str(1-chi2cdf(Dev,24-numParamsFuller),'%4.4f')])
        end
    end

    question = sprintf('\nDo you wish to compare the lesser and fuller model by way\nof the Likelihood Ratio Test [y/n]? ');
    wish = input(question,'s');
    if strcmpi(wish,'y')
        B = input('Type desired number of simulations B (try low, things take a while): ');
        [TLR pTLR paramsL paramsF TLRSim converged] = ...
            PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, ...
            paramsLesser, B, PF, 'searchOptions',options,'fullerlapserates',lapseratesfuller,'fullerthresholds',...
            thresholdsfuller,'fullerslopes',slopesfuller,'fullerguessrates',guessratesfuller,'lesserlapserates',...
            lapserateslesser,'lesserthresholds',thresholdslesser,'lesserslopes',slopeslesser,'lesserguessrates',...
            guessrateslesser,'lapseLimits',[0 1],'lapseFit',lapseFit);

        disp([sprintf('\n') 'Model comparison p-value by Monte Carlo: ' num2str(pTLR,'%4.4f')])
        if exist('chi2pdf.m') == 2
            disp(['Model comparison p-value by chi-square approximation: ' num2str(1-chi2cdf(TLR,numParamsFuller-numParamsLesser),'%4.4f')])
        end
    end
end

%what are the output and where are the?  I guees run and see? But I don't
%think this generate independent lapse rate fits for the right and left
%sides together. This may be fixed by having them lapse 'unconstrained' in
%the models. 