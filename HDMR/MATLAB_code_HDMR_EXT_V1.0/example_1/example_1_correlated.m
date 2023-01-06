%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%          GSA-HDMR_EXT: Linear Function with Correlated Inputs         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clc; clear; close all hidden;
% clear workspace and figures

% How many parameters?
d = 4; 
% Number of samples to be used
N = 5000;
% % Specify mean of samples
mu = 0.5*ones(1,d);
% All conditions are the same except the covariance matrix is now Eq. 44
C = eye(d); 
C(1,2) = 0.6; C(2,1) = 0.6; 
C(1,3) = 0.2; C(3,1) = 0.2;
C(3,2) = 0.2; C(2,3) = 0.2; 
% x values have mean of 0.5 and Sigma C
X = mvnrnd(mu,C,N); 
% X = (X - min(X))./ (max(X) - min(X));
% Model
Y = sum(X,2); 

% Random error
sigma2 = var(Y)/100; 
% Add random error to the model output
Y = Y + normrnd(0,sqrt(sigma2),N,1); 

% Specify HDMR_EXT options
options = struct('graphics',0,'basis',1,'maxorder',2,'m',3,'K',100,'R',3000,...
    'alfa',0.99,'method',1,'tolopt',1e-3);
% Now run the HDMR_EXT toolbox
[S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);