% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        1111   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           11 11   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       11  11   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE              11   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          11   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

clc; clear; close all hidden;           % clear workspace and figures

d = 5;                                  % # parameters?
N = 5000;                               % # samples to be used
mu = 0.5 * ones(1,d);                   % Sample mean, µ
C = eye(d);                             % Covariance matrix,Σ
X = mvnrnd(mu,C,N);                     % Draw N samples from N(µ,Σ) 
y = sum(X(:,1:d),2);                    % y = f(x) function
sigma2 = var(y)/100;                    % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);    % Add random error, y = Nx1 vector

options = struct('graphics',1, ...      % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1','bf2',0,'bf3',0,'m',2,...
    'K',100,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[S,Ss,Fx,Em,Xy] = HDMR(X,y,options);    % Now run the HDMR toolbox

[y_e,Y_ebf] = HDMR_model(Em,Xy);        % Look at last emulator
    
% options = struct('graphics',0,'maxorder',3,'maxiter',100,...
%     'm',2,'K',10,'R',300,'method',1,'alfa',0.99, ...
%     'lambda',10,'vartol',1e-3,'refit',1);