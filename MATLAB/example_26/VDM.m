function u = VDM(t,u_in,par,parameterization)
% Coupled ordinary differential equations of vandermeer population model
% Written by Jasper A. Vrugt
% University of California Irvine
% 
% The following papers are relevant
%  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and J.A.
%      Vrugt (2018), Probing the limits of predictability: data 
%      assimilation of chaotic dynamics in complex food webs, Ecology
%      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876
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


% Assign the states
P1 = u_in(1); P2 = u_in(2); Z1 = u_in(3); Z2 = u_in(4);

switch parameterization
    case 1
        beta = par(1);                      % predator coefficient
        alfa1 = par(2); alfa2 = alfa1;      % prey coefficient
        r1 = par(3);                        % growth rate 1st prey
        r2 = par(4);                        % growth rate 2nd prey
        g = par(5);                         % consumption rate
        K1 = par(6); K2 = K1;               % carrying capacity
        m1 = par(7); m2 = m1;               % mortality rate
        H = par(8);                         % parameter functional response
    case 2
        beta = par(1);                      % predator coefficient
        alfa1 = par(2);                     % prey 1 coefficient
        alfa2 = par(3);                     % prey 2 coefficient
        r1 = par(4);                        % growth rate 1st prey
        r2 = par(5);                        % growth rate 2nd prey
        g = par(6);                         % consumption rate
        K1 = par(7);                        % carrying capacity 1
        K2 = par(8);                        % carrying capacity 2
        m1 = par(9);                        % mortality rate 1
        m2 = par(10);                       % mortality rate 2      
        H = par(11);                        % parameter functional response
    otherwise
        error('VDM:unknown parameterization')
end

% Evaluate each of the four ordinary differential equations
dP1dt = r1*P1*(1 - (1/K1) * (P1 + alfa1 * P2)) - ...
    ( g*P1*Z1/( H + ( P1 + beta * P2 ) ) + ...
    g * beta * P1 * Z2/( H + ( beta * P1 + P2 ) ) );
dP2dt = r2*P2*(1 - (1/K2) * (alfa2 * P1 + P2)) - ...
    ( g*beta*P2*Z1/( H + ( P1 + beta * P2 ) ) + ...
    g * P2 * Z2/( H + ( beta * P1 + P2 ) ) );
dZ1dt = g * ( P1 + beta * P2 ) * Z1 / ...
    ( H + (P1 + beta*P2) ) - m1 * Z1;
dZ2dt = g * ( beta * P1 + P2 ) * Z2 / ...
    ( H + (beta*P1 + P2) ) - m2 * Z2;

% Assemble the dP/dt and dZ/dt in one vector u
u = [ dP1dt ; dP2dt ; dZ1dt ; dZ2dt ];

end
