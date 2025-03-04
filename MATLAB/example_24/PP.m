function u = PP(t,u_in,par)
% Coupled ordinary differential eqs. 1-predator-1-prey population model
% Written by Jasper A. Vrugt
% University of California Irvine

x = u_in(1); y = u_in(2);   % Assign state variables

r = par(1);         % intrinsic rate of prey natural increase
alfa = par(2);      % prop. constant linking prey mortality 
                    % to # prey and predators
m = par(3);         % mortality rate of predators
theta = par(4);     % prop. constant linking increase in predators 
                    % to # predators and prey
K = 50;             % maximum # preys environment can support 
                    % carrying capacity?

dxdt = r*x*(1 - x/K) - alfa*x*y;    % 1st of 2 ordinary differential eqs.
dydt = -m*y + theta*x*y;            % 2nd of 2 ordinary differential eqs.
u = [ dxdt ; dydt ];                % Assemble dx/dt & dy/dt in vector u

end
