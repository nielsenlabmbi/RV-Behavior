% 9/6/2023 - Modified from the PAL_PFBA_Demo in the Palamedes packages
%Aim is to fit curves onto free moving ferret behavior. Basically just a
%heavily edited demo as a starting point.

%RV "I am not very good at this."
%PAL_PFBA_Demo  Demonstrates use of Palamedes routine PAL_PFBA_Fit to fit a
%Psychometric Function (PF) to some data using a Bayesian criterion and
%determine the standard error estimates of the PF's parameters.

%
%Demonstrates usage of Palamedes function:
%-PAL_PFBA_Fit
%secondary:
%PAL_Logistic
%PAL_pdfNormal
%PAL_CumulativeNormal
%

%The inputs are:
% a logistic function
% @PAL_Logistic, @PAL_CumulativeNormal or other
PF = @PAL_CumulativeNormal ;

%Inputs from dotStairs.m
%Stimulus intensities
%StimLevels

%Number of trials at each entry of 'StimLevels'
%OutOfNum

%Number of positive (Right sided) responses
%   entries of 'StimLevels'  
%NumPos 

%Define which parameter values to include in prior/posterior. This 
%particular net is wide. As a rule, posterior should be (effectively) 
%contained within the grid, otherwise limits on grid will noticeably affect 
%estimates.
%orginal alpha and beta in demo
% grid.alpha = linspace(0,.12,201); %Threshold values, where can they fit between
% grid.beta = linspace(0,3,181);  % possible slopes, log-transformed values for beta

%alpha and beta from demo in PAL_PFBA_Fit.m
grid.alpha = [-180:5:180];  %if alpha is threshold then should this be the direction range?
grid.beta = [-.5:.01:.5];
grid.gamma = 0.5;               %using fixed value for guess rate ...
grid.lambda = 0.02;             %... and lapse rate, can be fixed like this or free, can be differnt for each side which can come from fitting

%Define a prior distribution across parameter space (optional). Seems
%appropriate for working in baysiean space
[a b g l] = ndgrid(grid.alpha,grid.beta,grid.gamma,grid.lambda);
prior = PAL_pdfNormal(a,0.06,0.03).*PAL_pdfNormal(b,1.5,1);
prior = prior./sum(sum(sum(sum(prior))));

%Fit function  - 9/6/23 Was having errors of 'Too many input arguments'
%with prior in the input list so I removed it for now. It has something to
%do with cutting out the varargin in PAL_PFBA_Fit when that was causing
%issues
[paramsValues2D posterior2D] = PAL_PFBA_Fit(StimLevels, NumPos, OutOfNum,grid, PF);
%got it to run with varargin but now have error "Output argument
%"paramsValues" (and maybe others) not assigned during call to
%"PAL_PFBA_Fit"."  where they had been beofore when running w/o prior
%arguments.  

%Fit again using free lapse rate
grid.lambda = [0:.001:.06];

% %Define a prior distribution across parameter space (optional)
[a b g l] = ndgrid(grid.alpha,grid.beta,grid.gamma,grid.lambda);
prior = PAL_pdfNormal(a,0.06,0.03).*PAL_pdfNormal(b,1.5,1).*l.^2.*(1-l).^98; %last two terms define beta distribution (minus normalization) with mode 0.02 on lapse rate
prior = prior./sum(sum(sum(sum(prior))));   %normalization happens here

%Fit function (RV also removed 'prior', prior from this list or arguments
[paramsValues3D posterior3D] = PAL_PFBA_Fit(StimLevels, NumPos, OutOfNum,grid, PF);

%Put summary of results w/ fixed lapse rate on screen
message = sprintf('\rLapse fixed at 0.02:');
disp(message);
message = sprintf('\rThreshold estimate: %6.4f',paramsValues2D(1,1));
disp(message);
message = sprintf('log10(Slope) estimate: %6.4f',paramsValues2D(1,2));
disp(message);
message = sprintf('Standard error of Threshold: %6.4f',paramsValues2D(2,1));
disp(message);
message = sprintf('Standard error of log10(Slope): %6.4f',paramsValues2D(2,2));
disp(message);

%Put summary of results w/ free lapse rate on screen
message = sprintf('\rLapse free:');
disp(message);
message = sprintf('\rThreshold estimate: %6.4f',paramsValues3D(1,1));
disp(message);
message = sprintf('log10(Slope) estimate: %6.4f',paramsValues3D(1,2));
disp(message);
message = sprintf('Lapse rate estimate: %6.4f',paramsValues3D(1,4));
disp(message);
message = sprintf('Standard error of Threshold: %6.4f',paramsValues3D(2,1));
disp(message);
message = sprintf('Standard error of log10(Slope): %6.4f',paramsValues3D(2,2));
disp(message);
message = sprintf('Standard error of lapse rate: %6.4f',paramsValues3D(2,4));
disp(message);


%% fixed lapse rate
%Figure showing posterior and fit w/ fixed lapse rate
f1 = figure('name','Bayesian Psychometric Function Fitting (fixed lapse rate)','units','pixels','position',[100 100 800 500]);
%colormap(f1,'jet')

%Contour plot of full posterior
ax1 = axes;
set(gca, 'units','pixels','position',[50 75 375 375]);
contour(gca,posterior2D')
set(gca, 'Xtick',[1:50:201],'FontSize',12);
set(gca, 'XtickLabel', {'0','.03','.06','.09','.12'});
xlabel(gca,'alpha');
set(gca, 'Ytick',[1:60:181],'FontSize',12);
set(gca, 'YtickLabel', {'0','1','2','3'});
ylabel(gca,'Log10(beta)');

%Data with fitted function
ax2 = axes;
hold on;
plot(gca, 0:.001:.12,PF([paramsValues2D(1,1) 10.^paramsValues2D(1,2) paramsValues2D(1,3) paramsValues2D(1,4)],0:.001:.12),'-','color',[0 .7 0],'linewidth',2) 
plot(gca, StimLevels,NumPos./OutOfNum,'ko','markersize',10,'markerfacecolor','k');
set(gca, 'units','pixels','position',[500 75 275 150]);
set(gca,'xlim',[0 .12],'ylim',[0.5 1],'xtick',StimLevels,'ytick',[0.5 1.0])
xlabel(gca,'Stimulus Intensity')
ylabel(gca,'proportion correct')

%Marginal posterior across threshold
ax3 = axes;
hold on;
plot(gca, grid.alpha,sum(posterior2D,2),'-') 
set(gca, 'units','pixels','position',[500 390 275 60]);
set(gca,'xlim',[0 .12],'ylim',[0 1.2*max(sum(posterior2D,2))],'xtick',StimLevels,'ytick',[])
xlabel(gca,'alpha')

%Marginal posterior across slope
ax4 = axes;
hold on;
plot(gca, grid.beta,sum(posterior2D,1),'-') 
set(gca, 'units','pixels','position',[500 280 275 60]);
set(gca,'xlim',[0 3],'ylim',[0 1.2*max(sum(posterior2D,1))],'xtick',[0:1:3],'ytick',[])
xlabel(gca,'log10(beta)')

%% Freelapse rate
%Figure showing posterior and fit w/ free lapse rate
f2 = figure('name','Bayesian Psychometric Function Fitting (free lapse rate)','units','pixels','position',[100 100 1000 600]);
%colormap(f2,'jet');

posterior3D = squeeze(posterior3D);

ax5 = axes;
[a, b, l] = ndgrid(grid.alpha, grid.beta, grid.lambda);
a = permute(a, [2 1 3]);
b = permute(b, [2 1 3]);
l = permute(l, [2 1 3]);
slice(gca,a,b,l,PAL_Scale0to1(permute(posterior3D,[2 1 3]))*64,[],[],[0:.01:.06]);
set(gca, 'units','pixels','position',[75 75 375 475]);
shading flat;
if ~exist('OCTAVE_VERSION');
    alpha('color');
    am = linspace(.25,1,64);
    alphamap(am);
end
axis(gca,[0 .12 0 3 0 0.06])
set(gca,'xtick',[.01:.02:.11]);
set(gca,'ytick',[0:1:3],'yticklabel',{'0','1','2','3'});
set(gca,'xgrid', 'off','ygrid','off','zgrid','off')
xlabel(gca,'alpha','fontsize',12)
ylabel(gca,'log10(beta)','fontsize',12)
zlabel(gca,'lambda','fontsize',12)

ax6 = axes;
image(grid.alpha, grid.beta,PAL_Scale0to1(flipud(sum(posterior3D,3)'))*64);
set(gca, 'units','pixels','position',[500 450 125 125]);
set(gca,'xtick',[.01:.05:.11]);
set(gca,'ytick',[0:1:3],'yticklabel',{'3','2','1','0'});
xlabel(gca,'alpha');
ylabel(gca,'log10(beta)');

ax7 = axes;
image(grid.alpha, grid.lambda,PAL_Scale0to1(flipud(squeeze(sum(posterior3D,2))'))*64);
set(gca, 'units','pixels','position',[675 450 125 125]);
set(gca,'xtick',[.01:.05:.11]);
set(gca,'ytick',[0:.02:.06],'yticklabel',{'.06','.04','.02','0'});
xlabel(gca,'alpha');
ylabel(gca,'lambda');

ax8 = axes;
image(grid.beta, grid.lambda,PAL_Scale0to1(flipud(squeeze(sum(posterior3D,1))'))*64);
set(gca, 'units','pixels','position',[850 450 125 125]);
set(gca,'xtick',[0:1:3]);
set(gca,'ytick',[0:.02:.06],'yticklabel',{'.06','.04','.02','0'});
xlabel(gca,'log10(beta)');
ylabel(gca,'lambda');

ax9 = axes;
plot(gca,grid.alpha,sum(sum(posterior3D,2),3));
set(gca, 'units','pixels','position',[500 325 125 75]);
set(gca,'xlim',[0 .12],'ylim',[0 1.2*max(sum(sum(posterior3D,2),3))])
set(gca,'xtick',[.01:.05:.11],'ytick',[]);
xlabel(gca,'alpha');

ax10 = axes;
plot(gca,grid.beta,sum(sum(posterior3D,1),3));
set(gca, 'units','pixels','position',[675 325 125 75]);
set(gca,'xlim',[0 3],'ylim',[0 1.2*max(sum(sum(posterior3D,1),3))])
set(gca,'xtick',[0:1:3],'ytick',[]);
xlabel(gca,'log10(beta)');

ax11 = axes;
plot(gca,grid.lambda,squeeze(sum(sum(posterior3D,1),2)));
set(gca, 'units','pixels','position',[850 325 125 75]);
set(gca,'xlim',[0 .06],'ylim',[0 1.2*max(sum(sum(posterior3D,1),2))])
set(gca,'xtick',[0:.02:.06],'ytick',[]);
xlabel(gca,'lambda');

%Data with fitted function
ax12 = axes;
hold on;
plot(gca, 0:.001:.12,PF([paramsValues2D(1,1) 10.^paramsValues2D(1,2) paramsValues2D(1,3) paramsValues2D(1,4)],0:.001:.12),'-','color',[0 .7 0],'linewidth',2) 
plot(gca, StimLevels,NumPos./OutOfNum,'ko','markersize',10,'markerfacecolor','k');
set(gca, 'units','pixels','position',[550 50 400 200]);
set(gca,'xlim',[0 .12],'ylim',[0.5 1],'xtick',StimLevels,'ytick',[0.5 1.0])
xlabel(gca,'Stimulus Intensity')
ylabel(gca,'proportion correct')