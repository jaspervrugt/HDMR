%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%           GSA-HDMR: Ishigami Function with Uncorrelated Inputs        %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clc; clear; close all hidden;
% clear workspace and figures

% How many parameters? 
d = 3; 
% Number of samples to be used
N = 8000;
% Latin hypercube
X = Latin(zeros(1,d),ones(1,d),N); % Note: much faster than built-in lhsdesign function
% Model
Y = sin(2*pi*X(:,1) - pi) + 7*(sin(2*pi*X(:,2) - pi)).^2 + 0.1*(2*pi*X(:,3) - pi).^4 .* sin(2*pi*X(:,1) - pi);


% Specify HDMR options
options = struct('graphics',0,'maxorder',2,'maxiter',100,'bf1',1','bf2',0,'bf3',0,'m',5,...
    'K',100,'R',4000,'method',1,'alfa',0.99,'lambda',0.10,'vartol',1e-3,'refit',1);
% Now run the HDMR toolbox
[S,Ss,Fx,Em,XY] = HDMR(X,Y,options);