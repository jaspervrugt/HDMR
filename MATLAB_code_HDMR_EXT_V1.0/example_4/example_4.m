%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script conduct global sensitivity analysis on a specific %
% dynamic rainfall-runoff model, using High Dimensional Model   %
% Representation with extended bases (HDMRext)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all
addpath(genpath([pwd '\samples']))
% Name of the M-file 
model_ID = 2002;
x_MAP_str = ['x_MAP_',num2str(model_ID),'.mat'];
x_post_str = ['x_post_',num2str(model_ID),'.mat'];
Y_MAP_str = ['Y_MAP_',num2str(model_ID),'.mat'];
Y_par_str = ['Y_par_',num2str(model_ID),'.mat'];
Y_tot_str = ['Y_tot_',num2str(model_ID),'.mat'];

% Load the posterior samples
load(x_MAP_str); load(x_post_str); load(Y_MAP_str);
load(Y_par_str);  % load(Y_tot_str);

% Number of observations
N_obs = size(Y_MAP,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%    1. Conduct Global SA using the posterior samples      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Define HDMR parameters and options
N = 10000;                       % Number of samples chosen
max_order = 2;                  % Maximum order of the component functions
options = struct('graphics',0,'basis',1,'maxorder',max_order,'maxiter',100,'bf1',1,'m',3,...
    'K',1,'R',N,'method',1,'alfa',0.99,'lambda',0.10,'vartol',1e-3,'refit',1);

% Dimension of the rainfall-runoff model
d_err = 4;                      % Dimension of the error model
d = size(x_post,2) - d_err;     % Dimension of the rainfall-runoff model

% Pre-allocation for sensitivity indices
switch max_order
    case 2
        n_ind = d + nchoosek(d,2);
    case 3
        n_ind = d + nchoosek(d,2) + nchoosek(d,3);
end
S_a = zeros(N_obs,n_ind);       % Structural index
S_b = zeros(N_obs,n_ind);       % Correlative index
S_t = zeros(N_obs,d);           % Total index
S_sum = zeros(N_obs,1);

% Randomly select N samples from the x_post
samp = randsample(size(x_post,1),N);
X = x_post(samp,1:d);           % Input matrix, X
X(1,:) = x_MAP(:,1:7);
Y_em = nan(N_obs,1);

% Start the loop
wtb = waitbar(0,'HDMRext computing...');
for i = 1:N_obs
    Y = Y_par(samp,i); 
    % Now run the HDMR_EXT toolbox
    [S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);
    % Extract the sensitivity indices
    for k = 1:n_ind
        S_a(i,k) = str2double(cell2mat(S(k+1,3)));
        S_b(i,k) = str2double(cell2mat(S(k+1,5)));
    end
    for k = 1:d
        S_t(i,k) = str2double(cell2mat(S(k+1,9)));
    end
    S_sum(i,:) = str2double(cell2mat(S(30,7)));
    % Extract the emulator
    Y_em(i,:) = Em.Y_e(1,:);
    clc
    waitbar(i/N_obs,wtb)
end
delete(wtb)

% Store the sensitivity indices to a struct variable
S_indices = struct('St',S_t,'Sa',S_a,'Sb',S_b);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% 2. Plot the structural and correlative sensitivity %%
%%            indices as a function of time           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% load state variables
load State_Variable.mat

plot_State_S_time_few(S_indices,S_sim,Y_MAP)

%% What if we use uncorrelated posterior/prior samples?
% Following step 1, we can also derive the sensitivity indices
% from posterior/prior hypercube by

% 1. Delineate the posterior/prior hypercube to get design matrix, X

% 2. Function evaluation to get matrix of Y

% 3. Rerun step 1 using the new 'X' and 'Y' as HDMRext input 

% We have included the prior samples, (X_prior, Y_par_prior)
% and posterior hypercube samples, (X_post_unc, Y_par_unc)
% in the folder. Feel free to load and run HDMRext using these samples!

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% 3. Plot the main effect derived from prior hypercube, %%
%%     posterior hypercube, and posterior distribution   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Now, we directly load the estimated results of those sensitivity indices
load('S_pos_correlated.mat')
load('S_pos_uncorrelated.mat')
load('S_prior_uncorrelated.mat')

plot_S_prior_posterior(S_indices_prior, S_indices_unc, S_indices)







