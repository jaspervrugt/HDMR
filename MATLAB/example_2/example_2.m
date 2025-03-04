% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      222222   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          22 22    %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE         22     %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE           22      %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      222222   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

clc; clear; close all hidden;               % clear workspace and figures

d = 5;                                      % # parameters?
N = 50000;                                  % # samples to be used
mu = 0.5*ones(1,d);                         % Sample mean, µ
C = eye(d) + diag([0.6 0.2 0 0.2],-1) ...   % Covariance matrix,Σ: Eq. 44
        + diag([0.6 0.2 0 0.2],1) ... 
        + diag([0.2 0 0],2) ...
        + diag([0.2 0 0],-2);
% C = [ 1.0 0.6 0.2 0.0 0.0 ;
%       0.6 1.0 0.2 0.0 0.0 ;
%       0.2 0.2 1.0 0.0 0.0 ;
%       0.0 0.0 0.0 1.0 0.2 ;
%       0.0 0.0 0.0 0.2 1.0 ]'
X = mvnrnd(mu,C,N);                         % Draw N samples from N(µ,Σ) 
X = (X - repmat(min(X),N,1))./ ...          % Normalize X - values
    repmat(max(X)-min(X),N,1);              
y = sum(X,2);                               % y = f(x) function
sigma2 = var(y)/100;                        % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',1, ...          % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1','bf2',0,'bf3',0,'m',3,...
    'K',100,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[S,Ss,Fx,Em,Xy] = HDMR(X,y,options);        % Now run the HDMR toolbox

% options = struct('graphics',1,'maxorder',3,'maxiter',100, ...
%     'm',2,'K',10,'R',300,'method',2,'alfa',0.99, ...
%     'lambda',10,'vartol',1e-3,'refit',1);
