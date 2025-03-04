function [X,Y] = evaluate_PP(X,n,d,u0,ode_options,dt,t_max)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function evaluates the 1-predator-1-prey model                     %
%                                                                         %
%  SYNOPSIS                                                               %
%   [X,Y] = evaluate_PP(X,n,d,u0,ode_options,dt,t_max)                    %
%  where                                                                  %
%    X           Nxd matrix of parameter vectors                          %
%    n           lenght of simulation                                     %
%                                                                         %
% Please check the following papers                                       %
%  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and     % 
%      J.A. Vrugt (2018), Probing the limits of predictability: data      %
%      assimilation of chaotic dynamics in complex food webs, Ecology     %
%      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876      %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

N = size(X,1);          % # parameter vectors
K = 2;                  % # species

% initialize output - N trials with n simulated counts and K = 2 species
Y = nan(N,n,K); prev = 0;
% now do the N simulations in parallel or not
for i = 1:N
    % Built-in ode45 function solves 2 diff. equations from t=0 to t=335
    [~,u] = ode45(@(t,u) PP(t,u,X(i,1:d)),0:dt:t_max,u0,ode_options);
    % now store each species in their own Y matrix
    Y(i,:,:) = reshape(u,1,n,K);
    % print progress
    fprintf(1, repmat('\b',1,prev)); 
    % print line
    prev = fprintf('MODEL EVALUATIONS: %3.2f%% DONE',100*(i/N));
end

end
