function [ci_thr, ci_s, ci_data, data_boot, t_boot, sl_boot] = bootstrap_pf(StimLevels,NumPos, OutOfNum, N, searchGrid, paramsFree,PF,options)
% BOOTSTRAP THRESHOLD AND CONFIDENCE INTERVALS
% number of simulations
fc = NumPos./OutOfNum;
boot_data = nan(length(NumPos),N);
data_boot = nan(length(NumPos),N);
for i = 1:length(StimLevels)
    tic
   %bootstrap the fraction correct
   boot_data(i,:) = binornd(OutOfNum(i),fc(i),[1 N]);
   data_boot(i,:) = boot_data(i,:)/OutOfNum(i);
end
t_boot = nan(1,N);
s_boot = nan(1,N);
sl_boot = nan(1,N);
for j = 1:N
    %perform fit for bootstrapped data
    if exist('options','var')
    paramsValues = PAL_PFML_Fit(StimLevels,boot_data(:,j)', ...
        OutOfNum,searchGrid,paramsFree,PF,'searchOptions',options,...
        'LapseLimits',[0 1],'GuessLimits',[0 1]);
    else
        paramsValues = PAL_PFML_Fit(StimLevels,boot_data(:,j)', ...
        OutOfNum,searchGrid,paramsFree,PF,'searchOptions',options,'LapseLimits',[0 1],'GuessLimits',[0 1]);
    end
    
    upperbound = 1-paramsValues(4); %paramsValues(4) is lapse rate, thus 1-paramsValues(4) is the upperbound on performance.
    guessrate = paramsValues(3); %percentage at 0
    
    middlepoint = PAL_CumulativeNormal(paramsValues,guessrate+(upperbound-guessrate)/2,'Inverse');
    sixtyseven = PAL_CumulativeNormal(paramsValues,guessrate+2*(upperbound-guessrate)/3,'Inverse');
    
    t_boot(j) = sixtyseven - middlepoint;
    sl_boot(j) = paramsValues(2);
    s_boot(j) = 1/t_boot(j);
    %get threshold for bootstrapped data
end
   %sort data into a distribution
   t_boot = sort(real(t_boot));
   s_boot = sort(real(s_boot));
   sl_boot = sort(real(sl_boot));
   %get 68% confidence intervals
   ci_thr = t_boot(round(N*[.16 .5 .84]));
   ci_s = s_boot(round(N*[.16 .5 .84]));
   sorted = sort(boot_data');
   %get 95% confidence intervals
   lims = sorted(round(N*[.025 .5 .975]),:);
   ci_data = [lims(1,:)./OutOfNum;      lims(2,:)./OutOfNum;       lims(3,:)./OutOfNum];
   toc

