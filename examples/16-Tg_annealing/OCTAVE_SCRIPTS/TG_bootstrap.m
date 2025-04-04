function [meantg,withinvar,trialcentertgs,rejected]=TG_bootstrap(densities,temperatures,numsyntheticsets,rejection,showfigs,runinparallel)

%% TG_bootstrap

%Version Notes; customized version of 7_01 for the book chapter

% The point of this function is to demonstrate details of the hyperbola fit,
% rejection criterion, and to yield histograms of Tg values that estimate
% the within-uncertainty (or sensitivity) of the hyperbola fit to
% short-timescale dynamical noise in the data.  This information will all
% be output in the form of several figures.  

% The script should be straightforward to use.  The variable "filename" is
% the filepath to the densities.txt file that contains the data set you
% would like to analyze.  The variable hypangle sets the hyperbolic angle
% used in the rejection criterion.  Other parameters are described in the
% "User inputs" section.  See also the Tg paper.  

% I use the following convention for variables that act as switches; 0 is
% off, 1 is on.

% *** NOTES on user inputs ***
% ``densities'' is a vector of densities; must be same dimensions as
% temperatures

%``temperatures'' is a vector of temperatures corresponding to the
%densities variable

%``numsyntheticsets'' is the number of synthetic datasets to generate

%``rejection'' is a flag, defaulted to 0 (or disabled) that runs the
%rejection criterion.  If rejection=0, the output variable rejected will
%always be set to 0 (i.e. not rejected)

%``showfigs'' is a flag that, when set to 0, disables all plotting.  When
%set to 1, the function generates various plots to illustrate the hyperbola
%fit and rejection criterion

%``runinparallel''; If set to 1, will run synthetic dataset computations in parallel.  You will need to set up parallel
%processing and enable it separately; otherwise the program will use the default
%parallel cluster object. Note that showraw2 (see below) is disabled if
%parallel is enabled

% *** Outputs ***
%``meantg'' is the mean value of T_g computed by the synthetic dataset
%analysis

%``withinvar'' is the within-simulation variance associated with the
%synthetic analysis

%``tghist'' is a vector of tg values arising from each synthetic dataset

%``rejected'' is a flag which is set to 0 if the set is not rejected and 1
%if it is


% User inputs and initializations

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------- User inputs ----------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
 



% *** NOTES TO USER; these parameters are needed for stability of the
% fitting process
% Set tgguess to a rough estimate of where you think Tg is
% I also find that scaleparam should be set to a value that is roughly
% where Tg is (or perhaps slightly below).
tgguess=300;            %guess value for Tg in Kelvin                 
scaleparam=300;        %this rescales the temperatures so that the characteristic domain of 
                        %the hyperbola is in a box with unit sides.
                        %this is needed to ensure that the fit is agnostic
                        %to the actual temperature range, thereby
                        %facilitating fitting

                        
% User inputs not included in the function call; can be changed here                        
scalepower=1.25;               %if greater than zero, terms in the least squares fit will be scaled by (i.e. divided by) temperature to this power
                              %Setting this to anything besides zero performs a weighted least squares
                            
tempscale=1.25;          %what should the temperature scaling of the noise model be?  Note, this parameter is different from 
                        %scalepower in that this is the scaling used to
                        %generate synthetic noise, not weighted least
                        %squares.  Note that residuals will be divided by temperature raised to the power tempscale                    
                       
                       
Percent_converged=0.80;             %greater than zero, less than one
                                    %how converged do you require the data
                                    %to be to the hyperbola asympototes?
                               
%debugging
residscalefactor=1;              %allows you to dial down or up the noise residuals; 1 is default               
showraw2=0;             %set to 1 to see individual hyperbola fits for synthetic data
            

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%----------------------- End User inputs ---------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

if nargin == 3
    rejection = 0;
    showfigs=0;
    runinparallel=0;
elseif nargin == 4
    showfigs=0;
    runinparallel=0;
elseif nargin == 5
    runinparallel=0;
end

%Definition for a hyperbola in terms of its six degrees of freedom

hyperbolafun=@(xdat,t)hyperbolafun_2(xdat,t);
weightedhyperbolafun=@(xdat,t) hyperbolafun(xdat,t)./((t).^(scalepower));            %only used for weighted least squares fits

%Guess parameters for the (rescaled) hyperbola fit
xo=[tgguess/scaleparam,1,1,1,1];

%Lower and upper bounds on the hyperbola fit parameters
lb = [ 0; 0; -3; -3; -30];
ub = [7; 10; 10; 10; 30];



%Re-organizing data and initializing a few variables
numtemps=numel(temperatures);
testtemps=min(temperatures):max(temperatures);

% Initial Hyperbola fit 


opts = optimset('Display','off','TolX',1e-12,'TolFun',1e-12);
if scalepower == 0          %i.e., if weight factors are set to 1
    [outs,resnormsx,ressx,exitflagx] = lsqcurvefit(hyperbolafun,xo,temperatures/scaleparam,densities,lb,ub,opts);
end

if scalepower > 0           %i.e. if weight factors are not 1
    [outs,resnormsx,ressx,exitflagx] = lsqcurvefit(weightedhyperbolafun,xo,temperatures/scaleparam,densities./((temperatures./scaleparam).^scalepower),lb,ub,opts);
end

%Hyperbola fit parameters
savedouts=outs;         %need to save the optimum parameters for use later in the linear noise propagation calculation

%vector of residuals and the sum of errors squared
rvec=densities - hyperbolafun(outs,temperatures/scaleparam);
resnormsx=sum(rvec.^2);

disp(' ')
disp(' ')
disp(' ')
disp('Starting analysis for current dataset')
disp('Tg from initial fit (in K)')
disp(outs(1)*scaleparam)
disp('Normed residual (smaller is better)')
disp(resnormsx)

% Rejection criterion
%The next few lines will use the hyperbola fits to determine where the
%asymptotic regimes are relative to the data.  
rejected=0;
if rejection == 1
    disp(strcat('Convergence parameter set to:',num2str(100*Percent_converged),'%'))

    bnds=hyperbolaAngle_2(outs',Percent_converged);
    lowerbound=bnds(1)*scaleparam;
    upperbound=bnds(2)*scaleparam;
    
    if showfigs == 1
        figure(1)
        tempouts=min([0.9*lowerbound,1.1*lowerbound,min(temperatures)]):max([1.1*upperbound,max(temperatures)]);
        plot(temperatures,densities,'bx')
        hold on
        plot(tempouts,hyperbolafun(outs,tempouts/scaleparam),'r')
        xlabel('Temperatures (in K)')
        ylabel('Density (g/cm^3)')
        title('Best fit hyperbola using all simulated data')
        plot([lowerbound,upperbound],[1,1],'g','LineWidth',3)
        plot(outs(1)*scaleparam,outs(2),'bo')
        plot([80,outs(1)*scaleparam],[outs(2)-outs(3)*(80/scaleparam-outs(1)),outs(2)],'k')
        legend('Simulated data','Best-fit hyperbola','Non-asymptotic regimes','Hyperbola Center','Asymptotes')
        plot([outs(1)*scaleparam,500],[outs(2),outs(2)-(outs(3)+outs(4))*(500/scaleparam-outs(1))],'k')
        hold off
    end
    
    if lowerbound < min(temperatures)
        rejected=1;
    end
    if upperbound > max(temperatures)
        rejected=1;
    end
    if rejected ==1
        disp('Data does not satisfy requirements for extracing Tg')
    else
        disp('Dataset not rejected')
    end

end

resids=(densities - hyperbolafun(outs,temperatures/scaleparam))*residscalefactor;         %residuals



% Noise modeling
trialsigma=std(resids./(temperatures.^tempscale));
trialnoise=zeros(numsyntheticsets,numtemps);
noisemat=random('norm',0,1,numsyntheticsets,max(size(temperatures)));
for jjjj=1:numsyntheticsets
    noisemat(jjjj,:)=noisemat(jjjj,:)*trialsigma.*(temperatures.^tempscale)';    %note that trialsigma doesn't appear in the lines below because it is contained in lowerchol       
    trialnoise(jjjj,:)=noisemat(jjjj,:);            %saves synthetic residuals to be recalled later
end

syntheticdensities=hyperbolafun(outs,temperatures/scaleparam) + trialnoise(1,:)';

if showfigs == 1
    
    %Plotting and analysis outputs
    figure(2)
    plot(temperatures,densities,'bx')
    hold on
    plot(testtemps,hyperbolafun(savedouts,testtemps/scaleparam),'r')
    hold off
    xlabel('Temperatures (in K)')
    ylabel('Density (g/cm^3)')
    title('Best fit hyperbola using all simulated data')
    
    figure(3)
    plot(temperatures,resids,'x')
    xlabel('Temperatures (in K)')
    ylabel('Density residuals (g/cm^3)')
    title('Residuals from best fit hyperbola')
    
    figure(4)
    plot(temperatures,resids./(temperatures.^tempscale),'x')
    xlabel('Temperatures (in K)')
    ylabel('Scaled residuals in units of Density/Temperature^{3/2} [g/(K^{3/2} \times cm^3)]')
    title('Scaled residuals')

    figure(5)
    plot(temperatures,noisemat(1,:),'x')
    xlabel('Temperatures (in K)')
    ylabel('Synthetic residuals based on noise model (g/cm^3)')
    title({'Example of synthetic residuals of best fit hyperbola'})
    
    figure(6)
    plot(temperatures,syntheticdensities,'bx');
    xlabel('Temperatures (in K)')
    ylabel('Synthetic Density (g/cm^3)')
    title({'Example of Synthetic Density vs Temperature curve','(Continuous curve is same as in Fig. 1 and not fit to synthetic data)'})
    hold on
    plot(testtemps,hyperbolafun(outs,testtemps/scaleparam),'r')
    hold off
end


trialcentertgs=zeros(1,numsyntheticsets);

if runinparallel == 0
    disp('Starting synthetic dataset analysis: progress indicated in 10% increments.')
    for jjjj=1:numsyntheticsets
        if mod(10*jjjj/numsyntheticsets,1) == 0
            disp(strcat('Synthetic set:',num2str(jjjj)))
        end
        syntheticdensities=hyperbolafun(savedouts,temperatures/scaleparam) + trialnoise(jjjj,:)';
        if scalepower == 0
            outs = lsqcurvefit(hyperbolafun,xo,temperatures/scaleparam,syntheticdensities,lb,ub,opts);
        end
        if scalepower > 0
            outs = lsqcurvefit(weightedhyperbolafun,xo,temperatures/scaleparam,syntheticdensities./((temperatures./scaleparam).^scalepower),lb,ub,opts);
        end
        if showraw2 == 1
            hold off
            plot(temperatures,syntheticdensities,'x')
            hold on
            plot(testtemps,hyperbolafun(outs,testtemps/scaleparam))
            plot([ex,ex],[0.9,1.5],'--g')
            axis([100, 600, 0.9,1.5]);
            pause
        end
        trialcentertgs(jjjj)=outs(1)*scaleparam;
    end
elseif runinparallel == 1
    disp('Starting synthetic dataset analysis: progress not indicated for parallel runs.')
    parfor jjjj=1:numsyntheticsets
        syntheticdensities=hyperbolafun(savedouts,temperatures/scaleparam) + trialnoise(jjjj,:)';
        if scalepower == 0
            outs = lsqcurvefit(hyperbolafun,xo,temperatures/scaleparam,syntheticdensities,lb,ub,opts);
        end
        if scalepower > 0
            outs = lsqcurvefit(weightedhyperbolafun,xo,temperatures/scaleparam,syntheticdensities./((temperatures./scaleparam).^scalepower),lb,ub,opts);
        end
        trialcentertgs(jjjj)=outs(1)*scaleparam;
    end
end

if showfigs == 1
    figure(7)
    hist(trialcentertgs)
    xlabel('T_g (in K)')
    ylabel('Number of synthetic sets')
    title('T_g values from from hyperbola center (entirely synthetic data)')
    hold on
end

disp('')
disp('')
disp('Mean hyperbola center T_g')
disp(mean(trialcentertgs));
disp('')
disp('Standard deviation of hyperbola center data')
disp(std(trialcentertgs));

meantg=mean(trialcentertgs);
withinvar=var(trialcentertgs);




end
