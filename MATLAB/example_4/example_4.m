% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE        4444   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE           44 44   %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       44  44   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          444444   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE          44   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

clc; clear; close all hidden;               % clear workspace and figures

d = 3;                                      % # parameters? 
N = 5000;                                   % # samples to be used
X = LH_sampling(zeros(1,3),ones(1,3),N);    % Latin hypercube sampling
X = (X - repmat(min(X),N,1))./ ...          % Normalize X - values
    repmat(max(X)-min(X),N,1);              
y = sin(2*pi*X(:,1) - pi) + ...             % y = f(x), function 
        7*(sin(2*pi*X(:,2) - pi)).^2 ...    
        + 0.1*(2*pi*X(:,3) - pi).^4 ...
        .* sin(2*pi*X(:,1) - pi);
sigma2 = var(y)/55;                         % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',0, ...          % Specify HDMR options
    'maxorder',3,'maxiter',100,'bf1',1', ...
    'bf2',0,'bf3',0,'m',4,'K',100, ...
    'R',500,'method',1,'alfa',0.01, ...
    'lambda',0.10,'vartol',1e-3,'refit',1);
[S,Ss,Fx,Em,Xy] = HDMR(X,y,options);        % Now run the HDMR toolbox
