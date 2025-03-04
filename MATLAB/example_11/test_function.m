function out = test_function(X)
% Modified Sobol test function
%
% REFERENCE
% Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto,
% M., and Tarantola, S. (2010). Variance based sensitivity
% analysis of model output. design and estimator for the
% total sensitivity index. Computer Physics Communications,
% 181, 259270.
%
% CONFIGURATION USED
% Campolongo, F., Saltelli, A., and Cariboni, J. (2011).
% From screening to quantitative sensitivity analysis. a
% unified approach. Computer Physics Communications,
% 182, 978988.

a = [100 0 100 100 100 100 1 10 0 0 9 0 100 100 4 100 100 7 100 2];
alpha = [1 4 1 1 1 1 0.4 3 0.8 0.7 2 1.3 1 1 0.3 1 1 1.5 1 0.6];
delt = [0.2942 0.2560 0.3004 0.5150 0.7723 0.4567 0.8390 0.1369 ...
    0.1558 0.4356 0.0257 0.3248 0.0718 0.9155 0.6877 0.5548 0.5835 ...
    0.8083 0.6309 0.8071];

for i = 1 : 20
    y(i) = ((1+alpha(i))*abs(2*(X(i)+delt(i)-fix(X(i)+delt(i)))-1) ...
        ^alpha(i)+a(i))/(1+a(i));
end

out = prod(y);