% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%   EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE      666666   %
%   EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE          66       %
%   EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE       666666   %
%   EE       XXXX   AAAAAA  MM   MM  PP      LL      EE          66  66   %
%   EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE      666666   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example from the following paper
%  Li, G., and H. Rabitz (2012), General formulation of HDMR component 
%      functions with independent and correlated variables, J. 
%      Math. Chem., 50, pp. 99-130

clc; clear; close all hidden;               % clear workspace and figures

d = 3;                                      % # parameters? 
N = 5000;                                   % # samples to be used
mu = [0.5 0.5 0.5];                         % Sample means, µ
s = [0.2 0.2 0.18];                         % Sample stds, σ
Cov = @(s,r12) [ s(1)^2 r12*s(1)*s(2) 0 ; ...
        r12*s(1)*s(2) s(2)^2 0 ; 0 0 s(3)^2 ];
r12 = 0.6; C = Cov(s,r12);                  % Create dxd covariance matrix, Σ
X = mvnrnd(mu,C,N);                         % Draw N samples from N(µ,Σ) 
X = (X - repmat(min(X),N,1))./ ...          % Normalize X - values
    repmat(max(X)-min(X),N,1);
a=[1 2]; b=[2 3]; c=[3 1 2]; e=[1 2 2 3];   % Inputs of user-defined function
y = user_function(X,mu,a,b,c,e);            % Evaluate function y = f(X)
sigma2 = var(y)/100;                        % Variance of random error
y = y + normrnd(0,sqrt(sigma2),N,1);        % Add random error, y = Nx1 vector

options = struct('graphics',1, ...          % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1,'bf2',0,'bf3',0,...
    'm',2,'K',10,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[S,Ss,Fx,Em,Xy] = HDMR(X,y,options);        % Now run the HDMR toolbox

% options = struct('graphics',1,'maxorder',3,'maxiter',100, ...
%     'm',2,'K',10,'R',300,'method',2,'alfa',0.99, ...
%     'lambda',10,'vartol',1e-3,'refit',1);

% subplot(2,2,1),plot(HDMR.X_n(:,1),Y_ebf(:,1),'b.');
% xlabel('x1','fontsize',16); ylabel('f1(x1)','fontsize',16); 
% set(gca,'fontsize',16);
% subplot(2,2,2),plot(HDMR.X_n(:,2),Y_ebf(:,2),'b.');
% xlabel('x2','fontsize',16); ylabel('f2(x2)','fontsize',16); 
% set(gca,'fontsize',16);
% subplot(2,2,3),plot(HDMR.X_n(:,3),Y_ebf(:,3),'b.');
% xlabel('x3','fontsize',16); ylabel('f3(x3)','fontsize',16); 
% set(gca,'fontsize',16);
% subplot(2,2,4),plot(X(:,3),Y_ebf(:,3),'b');
