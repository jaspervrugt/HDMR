function [y_res,y_ij,C2,T2,it] = HDMR_2nd(B2,y_res,R,n2,m2,vartol, ...
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
% Second order - individual estimation and backfitting                    %
%                                                                         %
%  SYNOPSIS                                                               %
%   [y_res,y_ij,C2,T2,it] = HDMR_2nd(B2,y_res,R,n2,m2,vartol, ...         %
%       lambda,maxiter,bf)                                                %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

C2 = zeros(m2,n2);          % Initialize coefficients
y_ij = zeros(R,n2);         % Initialize 2nd order contributions
T2 = zeros(m2,R,n2);        % Initialize T(emporary) matrix - 2nd
it = 0;                     % Initialize iteration counter

% Second order individual estimation
for j = 1:n2
    % Regularized least squares inversion ( minimize || C2 ||_2 )
    % C2(:,k) = ( B2(1:R,:,k)'*B2(1:R,:,k) + ...
    %               lambda * eye(m2) ) \ B2(1:R,:,k)' * y_res;
    B22 = B2(1:R,1:m2,j)'*B2(1:R,1:m2,j);
    if det(B22) == 0
        % Ill-defined --> C2(:,k) = 0;
    else
        T2(1:m2,1:R,j) = ( B22 + lambda * eye(m2) ) \ B2(1:R,1:m2,j)';
    end
    C2(1:m2,j) = T2(1:m2,1:R,j) * y_res;
    y_ij(1:R,j) = B2(1:R,1:m2,j) * C2(1:m2,j);
end

if (bf == 1)    % Second order backfitting
    var2b_old(1:n2) = sum(C2.^2,1); varmax = 1;
    while ( varmax > vartol ) && ( it < maxiter )
        for j = 1:n2
            y_r = y_res;
            for z = 1:n2
                if j ~= z
                    y_r = y_r - B2(1:R,1:m2,z) * C2(1:m2,z);
                end
            end
            % Regularized least squares inversion ( minimize || C2 ||_2 )
            % C2(:,k) = ( B2(:,:,k)'*B2(:,:,k) + ...
            %               lambda * eye(m2) ) \ B2(:,:,k)' * y_res;
            C2(1:m2,j) = T2(1:m2,1:R,j) * y_r;
        end
        var2b_new(1:n2) = sum(C2.^2,1);
        varmax = max(abs(var2b_new - var2b_old));
        var2b_old = var2b_new;
        it = it + 1;
    end
    % Now compute 2nd-order terms
    for j = 1:n2
        y_ij(1:R,j) = B2(1:R,1:m2,j) * C2(1:m2,j);
        % Subtract each first order term from residuals
        y_res = y_res - y_ij(1:R,j);
    end
end