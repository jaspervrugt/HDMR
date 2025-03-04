function [y_res,y_i,C1,T1,it] = HDMR_1st(B1,y_res,R,n1,m1,vartol, ...
    lambda,maxiter,bf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%          HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR              %
%          HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR             %
%          HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR            %
%          HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR            %
%          HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR             %
%          HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR              %
%          HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR             %
%          HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR             %
%          HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR            %
%          HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR           %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% First order - individual estimation and backfitting                     %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Y_res,Y_i,C1,T1,it] = HDMR_1st(B1,Y_res,R,n1,m1,vartol,...           %
%       lambda,maxiter,bf)                                                %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

C1 = zeros(m1,n1);      % Initialize coefficients
y_i = zeros(R,n1);      % Initialize 1st order contributions
T1 = zeros(m1,R,n1);    % Initialize T(emporary) matrix - 1st
it = 0;                 % Initialize iteration counter

% First order individual estimation
for j = 1:n1
    % Regularized least squares inversion ( minimize || C1 ||_2 )
    B11 = B1(1:R,1:m1,j)'*B1(1:R,1:m1,j);
    if det(B11) == 0
        % Ill-defined --> C1(1:m1,k) = 0;
    else
        T1(1:m1,1:R,j) = ( B11 + lambda * eye(m1) ) \ B1(1:R,1:m1,j)';
    end
    C1(1:m1,j) = T1(1:m1,1:R,j) * y_res;
    y_i(1:R,j) = B1(1:R,1:m1,j) * C1(1:m1,j);
end

if (bf == 1)    % First order backfitting
    var1b_old(1:n1) = sum(C1.^2,1); varmax = 1;
    while ( varmax > vartol ) && ( it < maxiter )
        for j = 1:n1
            y_r = y_res;
            for z = 1:n1
                if j ~= z
                    % Equivalent to statement (please test)
                    y_r = y_r - B1(1:R,1:m1,z) * C1(1:m1,z);
                end
            end
            % Regularized least squares inversion ( minimize || C1 ||_2 )
            % C1(:,k) = ( B1(:,:,k)'*B1(:,:,k) + ...
            %               lambda * eye(m1) ) \ B1(:,:,k)' * y_res;
            C1(1:m1,j) = T1(1:m1,1:R,j) * y_r;
        end
        var1b_new(1:n1) = sum(C1.^2,1);
        varmax = max(abs(var1b_new - var1b_old));
        var1b_old = var1b_new;
        it = it + 1;
    end
    % Now compute first-order terms
    for j = 1:n1
        y_i(1:R,j) = B1(1:R,1:m1,j) * C1(1:m1,j);
        % Subtract each first order term from residuals
        y_res = y_res - y_i(1:R,j);
    end
end

end
