%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%            GSA-HDMR: Linear Function with Uncorrelated Inputs         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clc; clear; close all hidden;
% clear workspace and figures

% How many parameters?
d = 4; 
% Number of samples to be used
N = 5000;
% % Specify mean of samples
mu = 0.5 * ones(1,d);
% All conditions are the same except the covariance matrix is now Eq. 44
C = eye(d);
% x values have mean of 0.5 and Sigma C
X = mvnrnd(mu,C,N);
% Model 
Y = sum(X(:,1:d),2);
% Random error
sigma2 = var(Y)/100;
% Add random error
Y = Y + normrnd(0,sqrt(sigma2),N,1);

% Specify HDMR options
options = struct('graphics',0,'maxorder',1,'maxiter',100,'bf1',1','bf2',0,'bf3',0,'m',2,...
    'K',100,'R',3000,'method',1,'alfa',0.99,'lambda',0.10,'vartol',1e-3,'refit',1);
% Now run the HDMR toolbox
[S,Ss,Fx,Em,XY] = HDMR(X,Y,options);