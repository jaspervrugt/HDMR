% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%       EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE           %
%       EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE               %
%       EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE            %
%       EE       XXXX   AAAAAA  MM   MM  PP      LL      EE               %
%       EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Please check the following papers
%  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and J.A.
%      Vrugt (2018), Probing the limits of predictability: data 
%      assimilation of chaotic dynamics in complex food webs, Ecology
%      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876
%  Vandermeer, J. (1993), Loose coupling of predator-prey cycles: 
%      entrainment, chaos, and intermittency in the classic MacArthur 
%      consumer-resource equations, American Naturalist, 141, pp. 687-716
%  Vandermeer, J. (2004), Coupled oscillations in food webs: balancing 
%      competition and mutualism in simple ecological models. American 
%      Naturalist, 163, pp. 857-867 
%  Vandermeer, J. (2006), Oscillating populations and biodiversity 
%      maintenance, Bioscience, 56, pp. 967-975

clc; clear; close all hidden;           % clear workspace and figures
model_version = 1;                      % Version of vandermeer population 
                                        % model (1: d=8 or 2: d=11 parameters)
N = 1000;                               % # VanDerMeer simulations? 
switch model_version                    % Vandermeer predator-prey model
    case 1  % d = 8 as in Massoud et al., EL, 2018
        d = 8;
        X_min = [.0001 .0001  0.01  0.01 0.1 0.50 0.001 1];
        X_max = [  1     2    2.50  2.50 2.5 2.50 0.700 3];
        label_par = {'$\beta$','$\alpha$','$r_{1}$','$r_{2}$', ...
            '$g$','$K$','$m$','$H$'};
    case 2  % d = 11 parameters (= test)
        d = 11;
        X_min = [ .0001 .0001 .0001 0.01  0.01 0.1 0.50 0.50 ...
            0.010 0.010 1 ];
        X_max = [   1     2     2   2.50  2.50 2.5 2.50 2.50 ...
            0.700 0.700 3 ];
        label_par = {'$\beta$','$\alpha_{1}$','$\alpha_{2}$', ...
            '$r_{1}$','$r_{2}$','$g$','$K_{1}$',...
            '$K_{2}$','$m_{1}$','$m_{2}$','$H$'};
end
X = LH_sampling(X_min,X_max,N);         % draw N parameter vectors
u0 = [ 1.0096 1.6658 0.5705 1.0006 ];   % initial count each four species
ode_options = odeset('abstol',1e-5, ... % ODE settings (trajectories subject to chaos!) 
    'reltol',1e-5, ...   
    'maxstep',0.1,'initialstep',0.1);
dt = 2; t_max = 1000; n = t_max/dt + 1; % Print step, start and end time

% Check whether file with model simulations exists or not
fid = fopen(['XY_',num2str(model_version),'.mat']);
% now check value of fid
switch fid
    case -1  % file does not exist
        % run model in parallel (requires distributed computing toolbox
        tbx = license('checkout','Distrib_Computing_Toolbox');
        % simulate count of each species at times (0,2,...,2000) days 
        [X,Y,SS] = VDM_model(X,n,d,u0,ode_options,dt,t_max,...
            model_version,tbx);
        % Also returns m = 3 summary stats simulations with n=501 counts.
        % This is in matrix SS(1:N,1:m,1:K), where K is number of species
        % Now save X and corresponding Y with subscript model version
        eval(strcat('save XY_',num2str(model_version),'.mat X Y SS'));
    otherwise
        % Load data instead
        eval(strcat('load XY_',num2str(model_version),'.mat'));     
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                NOW RUN HDMR CODE FOR EACH SIMULATION TIME
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

options = struct('graphics',1, ...      % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1,'bf2',0,'bf3',0,'m',5, ...
    'K',10,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[~,m,K] = size(SS); [~,n,~] = size(Y);  % # elements simulation and # species
for k = 1:K                             % Loop over each species    
    for m = 1:3                         % Loop over each summary metric
       [S,Ss,Fx,Em,XY] = ...            % Run HDMR toolbox 
           HDMR(X,Y,options);           % Need to downselect Y!
    end
end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% TO DO: try all kind of different ecologically relevant metrics 
% (from related literature). Then use HDMR to figure out which parameters 
% are most sensitive to what metric. This helps us to relate metrics to 
% parameters and see whether this makes sense in a time series
