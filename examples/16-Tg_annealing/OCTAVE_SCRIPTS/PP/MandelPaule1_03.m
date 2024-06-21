function [ darkvar,avg,varout ] = MandelPaule1_03(datain,ucin,lb,ub,pools,baseucin)

%Version 1_02
%Changes the tolerances in the nonlinear solver


%datain is a row vector of the input means
%ucin is a row vector of the variances associated with each mean
%pools is an optional matrix whose dimensions are as follows:
%   rows correspond to data pools, columns correspond to data sets in each
%   data pool.  The function will use the matrix pools to determine whether
%   there are any dependencies between the data pools, and adjust the
%   Mandel Paule algorithm to account for such dependencies.
%   if pools is invoked, then baseucin also has to be invoked.  This is the
%   set of uncertainties associated with the raw data sets (before
%   pooling).



numlabs=max(size(datain)) -1;
weight=@(s,y) (y + s).^(-1);
avgfun=@(y) sum(datain.*weight(ucin,y)) / sum(weight(ucin,y));
mpfun=@(y)sum(((datain - avgfun(y)).^2).*weight(ucin,y)) - numlabs;
varfun=@(y) ( sum(  ((datain - avgfun(y)).*weight(ucin,y)).^2) ) / (( sum(weight(ucin,y)) ).^2);


if nargin > 4
    dum=size(pools);
    numpools=dum(1);
    poolsize=dum(2);                       %size of each data pool
    dependencymatrix=zeros(numpools);       %a matrix that will determine the dependencies between the data pools
    pooluc=zeros(numpools);                %the cross uncertainties between data pools
    for pqrs=1:numpools
        for qwer=1:numpools
            intersection=intersect(pools(pqrs,:),pools(qwer,:));
            if min(size(intersection)) > 0
                dependencymatrix(qwer,pqrs)=max(size(intersection));
                pooluc(qwer,pqrs)=sum(baseucin(intersection));
            end
        end
    end
    pooluc=pooluc/(poolsize^2);
    dependencymatrix=dependencymatrix/poolsize;
    weight=@(s,y) (y + s).^(-1);
    avgfun=@(y) sum(datain.*weight(ucin,y)) / sum(weight(ucin,y));
    jj=1:numpools;
    genmpmat=@(y) weight(ucin(jj)',y)*weight(ucin(jj),y) .* (y .* dependencymatrix + pooluc); 
    genmpmodifier=@(y) sum(sum(genmpmat(y)));
    mpfun=@(y)sum(((datain - avgfun(y)).^2).*weight(ucin,y)) - numpools + genmpmodifier(y)/sum(weight(ucin,y));    
    varfun=@(y) genmpmodifier(y)/(sum(weight(ucin,y)).^2);
end

numtest=400;
testguess=zeros(1,numtest);
dumdumdum=zeros(1,numtest);
for uuuu=1:numtest
    testguess(uuuu)=mpfun(lb + (uuuu-1)*((ub-lb)/400));
    dumdumdum(uuuu)=lb + (uuuu-1)*((ub-lb)/400);
end

[dum,initguess]=min(abs(testguess));
darkvar=fzero(mpfun,(initguess-1)*((ub-lb)/400),optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8,'TolX',1e-10));
avg=avgfun(darkvar);
varout=varfun(darkvar);


end

