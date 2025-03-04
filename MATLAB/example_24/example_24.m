% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%       EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE           %
%       EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE               %
%       EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE            %
%       EE       XXXX   AAAAAA  MM   MM  PP      LL      EE               %
%       EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Example of the following paper
%  Cariboni, J., D. Gatelli, R. Liska, and A. Saltelli (2007), The role of 
%      sensitivity analysis in ecological modelling, Ecological Modelling, 
%      203, pp. 167–182.
%      http://izt.ciens.ucv.ve/ecologia/Archivos/ECO_POB%202007/ ...
%          ECOPO2_2007/Cariboni%20et%20al%202007.pdf
%  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and J.A.
%      Vrugt (2018), Probing the limits of predictability: data 
%      assimilation of chaotic dynamics in complex food webs, Ecology
%      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876

clc; clear; close all hidden;           % clear workspace and figures
d = 4;                                  % # model parameters
N = 100;                                % # pred-prey simulations? 

label_par = {'$r$','$\alpha$','$m$','$theta$'};
X_min = [ 0.8 0.2   0.6  0.05 ];        % min values of parameters    
X_max = [ 1.8 1.0   1.0  0.15 ];        % max values of parameters

X = LH_sampling(X_min,X_max,N);         % draw N parameter vectors
u0 = [ 8 2 ];                           % initial count each four species
ode_options = odeset('abstol',1e-5, ... % ODE settings (trajectories subject to chaos!) 
    'reltol',1e-5, ...   
    'maxstep',0.1,'initialstep',0.1);
dt = 1; t_max = 60; n = t_max/dt + 1;   % Print step, start and end time

% Check whether file with model simulations exists or not
fid = fopen('XY.mat','r');
% now check value of fid
switch fid
    case -1  % file does not exist
        % simulate count of each species at times (0,2,...,2000) days 
        [X,Y] = evaluate_PP(X,n,d,u0,ode_options,dt,t_max);
        % Now save X and corresponding Y with subscript model version
        save XY.mat X Y;
    otherwise
        % Load data instead
        load XY.mat X Y;     
end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                NOW RUN HDMR CODE FOR EACH SIMULATION TIME
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

options = struct('graphics',0, ...          % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1,'bf2',1,'bf3',1, ...
    'M',5,'K',10,'R',300, ...
    'method',1,'alfa',0.01, ...
    'lambda',0.10,'vartol',1e-3, ...
    'refit',1);
[N,n,K] = size(Y);                              
for k = 1:K                                 % Loop over each species
    for t = 1:10                            % Loop over time [limited to 10]
       [results(:,:,t,k),Ss,Fx, ...         % Run HDMR toolbox 
           Em,XY] = HDMR(X, ...             
           Y(1:N,t,k),options);             %#ok<SAGROW> 
       SS_interp{t,k} = Ss;                 %#ok
    end
end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% If you look at t = 9 and t = 61 then you get table 1 as in referenced 
% paper above. Results match quite closely. So for a paper we can start 
% with this case study and then followed by the two-predator-two-prey model
