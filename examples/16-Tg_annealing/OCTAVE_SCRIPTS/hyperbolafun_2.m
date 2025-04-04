function [ hx ] = hyperbolafun_2(hyperParam,x)
%HYPERBOLAFUN_2 Evaluates a parameterized hyperbola; see the chapter for details on this parametrization
%   
%   hx = hyperbolafun_2(hyperParam,x)
%

x0 = hyperParam(1);
y0 = hyperParam(2);
alpha = hyperParam(3);
gamma = hyperParam(4);
delta = hyperParam(5);

dx = x-x0;
hx = y0 - alpha*dx - gamma*(dx/2 + sqrt(dx.^2/4+exp(-delta)));

end


