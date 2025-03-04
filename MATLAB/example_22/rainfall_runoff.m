function [ Y ] = rainfall_runoff(x)
% Rainfall runoff model

% Store local variables in memory
persistent data y0 options tout

% Load the data and define local variables - only once
if isempty(tout)
    
    % Load the French Broad data
    daily_data = load('03451500.dly');
    % First two years are warm-up
    data.idx = [731:size(daily_data,1)]';
    % Define the PET and precipitation.
    data.P = daily_data(:,4); data.Ep = daily_data(:,5);
    %% Initial conditions
    y0 = 1e-5 * ones(5,1);
    %% Integration options
    options.InitialStep = 1;                 % initial time-step (d)
    options.MaxStep     = 1;                 % maximum time-step (d)
    options.MinStep     = 1e-5;              % minimum time-step (d)
    options.RelTol      = 1e-3;              % relative tolerance
    options.AbsTol      = 1e-3*ones(5,1);    % absolute tolerances (mm)
    options.Order       = 2;                 % 2nd order accurate method (Heun)
    %% Running time
    tout = [ 0 : size(data.P,1) ];
    
end

%% Assign parameters
data.Imax  = x(1);      % interception storage capacity (mm)
data.Sumax = x(2);      % unsaturated zone storage capacity (mm)
data.Qsmax = x(3);      % maximum percolation rate (mm/d)
data.aE    = x(4);      % evaporation coefficient
data.aF    = x(5);      % runoff coefficient
data.aS    = 1e-6;      % percolation coefficient
data.Kf    = x(6);      % fast-flow response time (d)
data.Ks    = x(7);      % slow-flow response time (d)

%% Run model C
y = crr_model(tout,y0,data,options);

% Now compute discharge
Y = diff(y(5,1:end))'; Y = Y(data.idx);

% Now calculate the summary metrics
% [ S_mod ] = calc_metrics( Y )';
% DB = S_mod;
% Now return in one vector Y and S_mod