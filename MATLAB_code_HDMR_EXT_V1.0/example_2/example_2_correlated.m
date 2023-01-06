%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%         GSA-HDMR_EXT: Ishigami Function with Correlated Inputs        %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clc; clear; %close all hidden;
% clear workspace and figures

% How many parameters? 
d = 3; 
% Number of samples to be used
N = 8000;
% Specify mean of samples
mu = 0.5*ones(1,d);
% Generate sample with correlation
pd=cell(1,3);
pd{1} = makedist('Uniform',0,1);
pd{2} = makedist('Uniform',0,1);
pd{3} = makedist('Uniform',0,1);
correlation = eye(3); correlation(1,3) = 0.5; correlation(3,1) = 0.5;
X = lhsgeneral(pd,correlation,N);

% Model
Y = sin(2*pi*X(:,1) - pi) + 7*(sin(2*pi*X(:,2) - pi)).^2 + 0.1*(2*pi*X(:,3) - pi).^4 .* sin(2*pi*X(:,1) - pi);

% Specify HDMR_EXT options
options = struct('graphics',0,'basis',1,'maxorder',2,'m',6,'K',100,'R',N/2,...
    'alfa',0.99,'method',1,'tolopt',1e-3);
% Now run the HDMR_EXT toolbox
[S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);