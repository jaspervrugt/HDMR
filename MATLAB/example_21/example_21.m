% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%       EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE           %
%       EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE               %
%       EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE            %
%       EE       XXXX   AAAAAA  MM   MM  PP      LL      EE               %
%       EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

clc; clear; close all hidden;       % clear workspace and figures

d = 4;                              % # parameters?
N = 500;                            % # samples to be used
solver = 'explicit';                % Numerical solution: 'explicit'/'implicit'

fid = fopen(['XY_',solver,'.mat']); % Load XY data
switch fid  % Does this file exist?
    case -1 % no --> sample X and create corresponding model output
        % Simulate the extent for temperature 500...700 for N parameter vectors
        [X,label_par,Y] = evaluate_TE(N,solver);
        % now save X and corresponding Y with subscript model version
        eval(strcat('save XY_',solver,'.mat X Y'));
    otherwise
        % XY data file exists so you can load it instead
        eval(strcat('load XY_',solver,'.mat'));
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%              NOW RUN HDMR_EXT CODE FOR EACH SIMULATION TIME
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

options = struct('graphics',0, ...      % Specify HDMR_EXT options - d = 2     
    'maxorder',2,'maxiter',100, ...     % --> maxorder cannot exceed 2
    'bf1',1,'bf2',0,'bf3',0,...
    'm',5,'K',50,'R',200,'method',1, ...
    'alfa',0.01,'lambda',0.10, ...
    'vartol',1e-3,'refit',1);
[N,n] = size(Y);                        % # simulations, # values per trial
                                        % Determine N, as "bad" simulations 
                                        % may arise due to numerical errors
results = cell(5,13,n);                 % initialize results
for t = 1:n                             % loop over each temperature
    [results(1:5,1:13,t),Ss,Fx, ...
        Em,Xy] = HDMR(X,Y(1:N,t),options);
end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% extract total sensitivity of x1,x2,x1x2 and sum
out = str2double(reshape(results(2:5,7,:),4,201));
% can now plot sensitivites as function of temperature
T = (500:1:700);
plot(T,out(1,:),'r','linewidth',2); hold on;
plot(T,out(2,:),'b','linewidth',2);
plot(T,out(3,:),'g','linewidth',2);
plot(T,out(4,:),'k','linewidth',2);
[objh] = legend({'$x_{1}$','$x_{2}$','$x_{1}/x_{2}$','$\sum$'}, ...
    'interpreter','latex','fontsize',18);
set(objh,'linewidth',3);
set(gca,'fontsize',16);
xlabel('Temperature','fontsize',18,'fontweight','bold');
ylabel('Sensitivity','fontsize',18,'fontweight','bold');
