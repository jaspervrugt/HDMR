function [y_ijk,C3,T3,it] = HDMR_3rd(B3,y_res,R,n3,m3,vartol,lambda, ...
    maxiter,bf)
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
% Third order - individual estimation and backfitting                     %
%                                                                         %
%  SYNOPSIS                                                               %
%   [y_ijk,C3,T3,it] = HDMR_3rd(B3,y_res,R,n3,m3,vartol,lambda, ...       %
%       maxiter,bf)                                                       %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

C3 = zeros(m3,n3);          % Initialize coefficients
y_ijk = zeros(R,n3);        % Initialize 3rd order contributions
T3 = zeros(m3,R,n3);        % Initialize T(emporary) matrix - 3rd
it = 0;                     % Initialize iteration counter

% Third order individual estimation
for j = 1:n3
    % Regularized least squares inversion ( minimize || C3 ||_2 )
    % C3(:,k) = ( B3(1:R,:,k)'*B3(1:R,:,k) + ...
    %               lambda * eye(m3 ) ) \ B3(1:R,:,k)' * y_res;
    B33 = B3(1:R,1:m3,j)'*B3(1:R,1:m3,j);
    if det(B33) == 0
        % Ill-defined --> C3(:,k) = 0;
    else
        T3(1:m3,1:R,j) = ( B33 + lambda * eye(m3) ) \ B3(1:R,1:m3,j)';
    end
    C3(1:m3,j) = T3(m3,1:R,j) * y_res;
    y_ijk(1:R,j) = B3(1:R,1:m3,j) * C3(1:m3,j);
end

if (bf == 1)    % Third order backfitting
    var3b_old(1:n3) = sum(C3.^2,1); varmax = 1;
    while ( varmax > vartol ) && ( it < maxiter )
        for j = 1:n3
            y_r = y_res;
            for z = 1:n3
                if j ~= z
                    y_r = y_r - B3(1:R,1:m3,z) * C3(1:m3,z);
                end
            end
            % Regularized least squares inversion ( minimize || C3 ||_2 )
            % C3(:,k) = ( B3(:,:,k)'*B3(:,:,k) + ...
            %               lambda * eye(m3) ) \ B3(:,:,k)' * y_res;
            C3(1:m3,j) = T3(1:m3,1:R,j) * y_r;
        end
        var3b_new(1:n3) = sum(C3.^2,1);
        varmax = max(abs(var3b_new - var3b_old));
        var3b_old = var3b_new;
        it = it + 1;
    end
    % Now compute 3rd-order terms
    for j = 1:n3
        y_ijk(1:R,j) = B3(1:R,1:m3,j) * C3(1:m3,j);
        % Subtract each first order term from residuals (= not needed)
        % y_res = y_res - yijk_bf(:,k);
    end
end