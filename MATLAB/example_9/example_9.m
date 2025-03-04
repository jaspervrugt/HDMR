% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      999999   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          99  99   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       999999   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              99   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      999999   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example 1 of the following paper
%  Chastaing, G., F. Gamboa, and C. Prieur (2012), Generalized Hoeffding-
%      Sobol decomposition for dependent variables - application to 
%      sensitivity analysis, Electron. J. Statistics, 6, pp. 2420-2448, 
%      https://doi.org/10.1214/12-EJS749

clc; clear; close all hidden;               % clear workspace and figures

d = 2;                                      % # parameters? 
N = 1000;                                   % # samples to be used
[mu1,mu2] = deal(zeros(d,1));               % sample means µ1 & µ2 normal mixture
w1 = 0.2;                                   % weight of first normal
r12 = 0.4;                                  % covariance X1,X2 2nd normal 
C1 = eye(d); C2 = [ 0.5, r12 ; r12, 0.5 ];  % covariance matrices, Σ1 & Σ2
X = user_draw(mu1,mu2,C1,C2,w1,N,d);        % Draw N samples X normal mixture
y = X(:,1) + X(:,2) + X(:,1) .* X(:,2);     % y = f(X) function
sigma2 = var(y)/100;                        % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',1, ...          % Specify HDMR options
    'maxorder',2,'maxiter',100, ...
    'bf1',1,'bf2',1,'bf3',1,...
    'm',2,'K',10,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[S,Ss,Fx,Em,Xy] = HDMR(X,y,options);        % Now run the HDMR toolbox
