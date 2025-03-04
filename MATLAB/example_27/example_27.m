% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%       EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE           %
%       EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE               %
%       EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE            %
%       EE       XXXX   AAAAAA  MM   MM  PP      LL      EE               %
%       EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clc; clear; close all hidden;               % clear workspace and figures

N = 1e5;                                    % Parameter vectors (= samples)
X = rand(N,3);                              % Sample from U(0,1) N times
f = @(x) x(:,1) + 2*x(:,2) ...              % y = f(x) function
        + 3*x(:,3) + x(:,1).*x(:,2);
y = f(X);                                   % Compute function output

options = struct('graphics',1, ...          % Specify HDMR options
    'basis',1,'maxorder',2,'m',3, ...
    'K',1,'R',N,'alfa',0.01, ...
    'method',1,'tolopt',1e-3);
[S,Ss,Fx,Em,Xy] = HDMR_EXT(X,y,options);    % Run HDMR toolbox
