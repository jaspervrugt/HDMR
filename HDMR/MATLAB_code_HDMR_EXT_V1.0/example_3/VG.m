function theta = VG(x,h)
theta_s = x(1); theta_r = x(2); alpha = x(3); n = x(4);
theta = theta_r + (theta_s - theta_r)*(1+(alpha*abs(h)).^n).^(1/n-1);