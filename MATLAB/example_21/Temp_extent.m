function [T_s,a_s,flag] = Temp_extent(x,plugin)
% Solves temperature - extent relationship for given E = x(1) and A = x(2)

flag = 1;       % "issues" flag
R = 8.314e-3;   % Universal gas constant
% Two numerical solution methods for stiff differential equation
switch lower(plugin.solution)
    case 'implicit'
        [T_s,a_s] = ode23s(@ode_model,plugin.T_e,...
            plugin.a0,plugin.options,x,plugin.beta,R);  % implicit solver
    case 'explicit' % eplicit solver -> fixed dTee
        nT = 10000;                                     % number of temperature steps      
        dT = (plugin.T_e(end) - plugin.T_e(1))/nT;      % integration step of temperature
        T = plugin.T_e(1):dT:plugin.T_e(end);         % simulated temperatures
        [T_s,a_s] = explicit_model(T,plugin.T_e,...
            plugin.a0,x,plugin.beta,R,dT,nT);           % explicit solver
    otherwise
        error('HDMR:Temp_extent:unknown integration')
end
% now check output - and set to zero if nan
if sum(isnan(a_s)) > 0, a_s = zeros(plugin.n,1); flag = -1; end

% -------------------------------------------------------------------------
% SECUNDAIRY FUNCTIONS: ode_model(T_e,a,x,beta,R)
% -------------------------------------------------------------------------
function da = ode_model(T_e,a,x,beta,R)
% implicit solver - return da; da = f(T,a,x)
da = (1/beta) * 10^x(2) * exp ( -x(1)/R./T_e ) .* (1 - a);
return

% -------------------------------------------------------------------------
% SECUNDAIRY FUNCTIONS: ecplicit_model(T_e,a0,x,beta,R,dT,nT)
% -------------------------------------------------------------------------
function [T_e,a_s] = explicit_model(T,T_e,a0,x,beta,R,dT,nT)
% explicit solver integrates with dT
a = nan(nT+1,1); a(1) = a0;                     % preallocate a and its initial state
for i = 2:nT+1
    da = (1/beta) * 10^x(2) * ...
        exp( -x(1)/R./T(i-1) ) .* (1 - a(i-1) );
    a(i) = a(i-1) + da * dT;
end
% Now interpolate [T,a] at T_e
a_s = interp1(T,a,T_e);
return