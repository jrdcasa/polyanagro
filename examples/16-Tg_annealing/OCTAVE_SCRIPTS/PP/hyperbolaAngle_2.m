function [ xBnds ] = hyperbolaAngle_2( hyperParam, r )
%The function hyperbolaAngle_2 Returns x at which the hyperbola slope is the fraction r and (1-r) of its two asymptotic values.
% The parameter 0<r<1 can be thought of as the percent convergence of the hyperbola to its asymptotes.
% Here we assume the hyperbola parametrization according to the function hyperbolafun_2.m
% hyperParam = [x0;y0;alpha;gamma;delta] then 
%
%   df(xLow) = m0 = -(alpha+gamma)+r*gamma
%   df(xHigh)= m1 = -(alpha+gamma)+(1-r)*gamma
%


x0 = hyperParam(1);		%the temperature coordinate of the hyperbola center
delta = hyperParam(5);	

dX   = exp(-delta/2)*(2*r-1)/sqrt(r*(1-r));
xLow = x0-dX;
xHigh= x0+dX;

xBnds = [xLow;xHigh];		%upper and lower bounds on what temperature are considered to be within the transition region, i.e. not those 					%temperatures for which the hyperbola is r-percent converged to its asymptotes.

end

