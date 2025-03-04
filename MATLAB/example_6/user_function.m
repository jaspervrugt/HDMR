function y = user_function(X,mu,a,b,c,e)
% Evaluate user defined function, y = f(X), for the Nxd matrix of
% parameter vectors X
% Written by Jasper A. Vrugt
% University of California Irvine

g1 = @(x,mu,a,b) (a(2)*(x(:,1) - mu(1)) ... 
        + a(1)).*(b(2)*(x(:,2) - mu(2)) + b(1));
g2 = @(x,mu,c) c(3)*(x(:,2) - mu(2)).^2 ...
        + c(2)*(x(:,2) - mu(2)) + c(1);
g3 = @(x,mu,e) e(4)*(x(:,3) - mu(3)).^3 ...
        + e(3)*(x(:,3) - mu(3)).^2 ...
        + e(2)*(x(:,3)-mu(3)) + e(1);
y = g1(X,mu,a,b) + g2(X,mu,c) + g3(X,mu,e);

end