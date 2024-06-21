function [ darkvar,avg,varout ] = Dark_MLE(datain,ucin,lb,ub)

%This function computes the weighted-mean average and dark uncertainty associated with estimates of an arbitrary quantity that come with within-simulation uncertainties.

%Inputs are as follows:
%1) datain is a row vector of inputs means of the quantity of interest
%2) ucin is a row vector of variances associated with each mean
%3) lb and ub are lower and upper bounds respectively on what we expect the dark uncertainty to be (note that dark uncertainty has units of 			%variance in this function
%4) Outputs are: darkvar, i.e. the dark uncertainty, avg, i.e. the weighted mean statistic, and varout, i.e. the uncertainty (in variance 		%units) of the average.


weight=@(s,y) (y + s).^(-1);			%weight function associated with the dark uncertainty / weighted-mean statistic
avgfun=@(y) sum(datain.*weight(ucin,y)) / sum(weight(ucin,y));			%weighted-mean average as a function of the (as-of-yet-unknown) dark uncertainty
mlefun=@(y)1.*sum(((datain - avgfun(y)).^2).*(weight(ucin,y).^2)) - 1.*sum(weight(ucin,y));	%MLE function to solve for y
varfun=@(y) ( sum(  ((datain - avgfun(y)).*weight(ucin,y)).^2) ) / (( sum(weight(ucin,y)) ).^2);	%uncertainty in the mean

%Preliminary test to find a good initial guess for the dark uncertainty
numtest=400;
testguess=zeros(1,numtest);
dumdumdum=zeros(1,numtest);
for uuuu=1:numtest
    testguess(uuuu)=mlefun(lb + (uuuu-1)*((ub-lb)/numtest));
    dumdumdum(uuuu)=lb + (uuuu-1)*((ub-lb)/numtest);
end



[dum,initguess]=min(abs(testguess));	%initial guess for the dark uncertainty

%Computation to find dark uncertainty (nonlinear solve)
darkvar=fzero(mlefun,(initguess-1)*((ub-lb)/numtest),optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8,'TolX',1e-10));

%Weighted-mean statistic
avg=avgfun(darkvar);

%Uncertainty in the weighted-mean statistic
varout=varfun(darkvar);

end

