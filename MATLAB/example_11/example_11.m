% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE    11  11     %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE        11  11     %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE     11  11     %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        11  11     %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE    11  11     %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clc; clear; close all hidden;               % clear workspace and figures

% Modified Sobol test function
d = 20;                                     % # parameters? 
N = 5000;                                   % # samples to be used
X = rand(N,d);                              % Sample X-values from U(0,1)
y = nan(N,1);                               % Initialize y values
for i = 1:N
    y(i,1) = test_function(X(i,1:d));       % Evaluate modified Sobol func
end

options = struct('graphics',1, ...          % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1,'bf2',1,'bf3',1, ...
    'm',2,'K',10,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[S,Ss,Fx,Em,Xy] = HDMR(X,y,options);        % Now run the HDMR toolbox
