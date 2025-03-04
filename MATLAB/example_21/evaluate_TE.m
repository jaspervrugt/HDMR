function [X,label_par,Y] = evaluate_TE(N,solver)
% This function runs the Temperature - Extent model N different times

% specify temperature range, beta, initial state and numerical solution
plugin = struct('T_e',500:1:700,'beta',1/6,'a0',0, ...
    'solution',solver,'n',201,'options',odeset);
% number of model parameters
d = 2;                  
% minimum values of parameters
X_min = [ 1.0 1.00 ];   
% maximum values of parameters
X_max = [ 500 20.0 ];   
% labels for parameters
label_par = {'$E$','$A$'};

% sample randomly the parameter values within their ranges
X = Latin(X_min,X_max,N);
% initialize model output - N trails with n simulated counts and k = 4 species
Y = nan(N,plugin.n); 
% set prev = 0 for printing
prev = 0; flag = nan(N,1);
% Now evaluate the model (model is fast so parfor not needed)
for i = 1:N
    [~,a_s,flag(i,1)] = Temp_extent(X(i,1:2),plugin);
    % now store each species in their own Y matrix
    Y(i,1:plugin.n) = a_s';
    % print progress
    fprintf(1, repmat('\b',1,prev)); prev = ...
        fprintf('MODEL EVALUATIONS: %3.2f%% DONE',100*(i/N)); 
end
% Now remove "bad" simulations
idx = (flag == 1);
Y = Y(idx,1:plugin.n);
% Same as with X
X = X(idx,1:d);