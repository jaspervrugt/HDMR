function [SSR] = VG_fit(x,h_meas,t_meas)
% Define function
ts = x(1); 
tr = x(2);
a = x(3);
n = x(4);
m = (n-1)/n;
% Compute simulated theta for given parameters and measured h-values
t_sim = tr + (ts - tr) * ( 1 + (a*abs(h_meas)).^n).^(-m);
% Now compute residuals ( = m x 1 vector ; m = 10 in our example)
res = t_meas - t_sim;
% Now compute sum of squared residuals
SSR = sum(res.^2);