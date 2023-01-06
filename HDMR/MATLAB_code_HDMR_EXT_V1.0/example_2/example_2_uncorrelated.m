%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%         GSA-HDMR_EXT: Ishigami Function with Uncorrelated Inputs      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clc; clear; close all hidden;
% clear workspace and figures

% How many parameters? 
d = 3; 
% Number of samples to be used
N = 8000;
% Latin hypercube
X = Latin(zeros(1,3),ones(1,3),N); % Note: much faster than built-in lhsdesign function
% Model
Y = sin(2*pi*X(:,1) - pi) + 7*(sin(2*pi*X(:,2) - pi)).^2 + 0.1*(2*pi*X(:,3) - pi).^4 .* sin(2*pi*X(:,1) - pi);

% Specify HDMR_EXT options
options = struct('graphics',0,'basis',1,'maxorder',2,'m',6,'K',100,'R',N/2,...
    'alfa',0.99,'method',1,'tolopt',1e-3);

% Now run the HDMR_EXT toolbox
[S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);