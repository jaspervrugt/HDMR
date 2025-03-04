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
%  Benincà, E., J. Huisman, R. Heerkloss, K.D. Jöhnk, P. Branco, E.H. Van
%      Nes, M. Scheffer, and S.P. Ellner (2008), Chaos in a long-term 
%      experiment with a plankton community, Nature, 451, pp. 822-826, 
%      doi:10.1038/nature06512, www.nature.com/doifinder/10.1038/nature06512
%  Benincà, E., K.D. Jöhnk, R. Heerkloss, and J. Huisman (2009), Coupled 
%      predator–prey oscillations in a chaotic food web, Ecology Letters, 
%      12, pp. 1367-1378, doi:10.1111/j.1461-0248.2009.01391.x
%  Vandermeer, J. (2006), Oscillating populations and biodiversity 
%      maintenance, Bioscience, 56, pp. 967-975
%  Vandermeer, J. (2004), Coupled oscillations in food webs: balancing 
%      competition and mutualism in simple ecological models, The American 
%      Naturalist, 163, pp. 857-867 
%  Vandermeer, J. (1993), Loose coupling of predator-prey cycles: 
%      entrainment, chaos, and intermittency in the classic MacArthur 
%      consumer-resource equations, The American Naturalist, 141, 
%      pp. 687-716

clc; clear; close all hidden;           % clear workspace and figures
model_version = 1;                      % Version of vandermeer population 
                                        % model (1: d=8 or 2: d=11 parameters)
N = 500;                                % # VanDerMeer simulations? 
switch model_version                    % Vandermeer predator-prey model
    case 1  % d = 8 as in Massoud et al., EL, 2018
        d = 8;
        X_min = [ 1e-4  1e-4  0.01  0.01 0.1 0.50 1e-3  1];
        X_max = [  1     2    2.50  2.50 2.5 2.50 0.700 3];
        label_par = {'$\beta$','$\alpha$','$r_{1}$','$r_{2}$', ...
            '$g$','$K$','$m$','$H$'};
    case 2  % d = 11 parameters (= test)
        d = 11;
        X_min = [ 1e-4  1e-4  1e-4  0.01  0.01 0.1 0.50 0.50 ...
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
dt = 3.35; t_max = 2656.55;             % Print step, start & end time
n = t_max/dt + 1;                       % # print steps
K = 4;                                  % # species

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
        % Now compute the residual with observed data
        XLS_sheet = readmatrix(['data_Beninca_et_al_Nature&' ...
            'Ecology_Letters.xlsx'],'sheet','sheet1','range','A1:E794');
        % Assign measured data
        MeasData = XLS_sheet(:,2:5)';
        % initialize matrix of residuals
        Res = nan(N,n,K);
        % Now loop to determine residuals
        for k = 1:K % Species
            for i = 1:N % Parameter vectors
                Res(i,1:n,k) = Y(i,1:n,k) - MeasData(k,1:n); 
            end
        end
        % Now save X and corresponding Y with subscript model version
        eval(strcat('save XY_',num2str(model_version),'.mat X Y SS Res'));
    otherwise
        % Load data instead
        eval(strcat('load XY_',num2str(model_version),'.mat'));     
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                NOW RUN HDMR CODE FOR EACH SIMULATION TIME
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

options = struct('graphics',0, ...      % Specify HDMR options
    'maxorder',3,'maxiter',100, ...
    'bf1',1,'bf2',0,'bf3',0,'m',5, ...
    'K',10,'R',300,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[~,n,K] = size(Res);                    % # elements simulation and # species
for k = 1:K                             % Loop over each species    
    for t = 2:10                        % Loop over each time [limited to 10]
       [results(:,:,t,k),Ss, ...        % Run HDMR toolbox
           Fx,Em,XY] = HDMR(X, ...
           Res(1:N,t,k),options);       % Need to downselect Y!
    end
end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% TO DO: try all kind of different ecologically relevant metrics 
% (from related literature). Then use HDMR to figure out which parameters 
% are most sensitive to what metric. This helps us to relate metrics to 
% parameters and see whether this makes sense in a time series