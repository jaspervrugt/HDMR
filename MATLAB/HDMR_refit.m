function [Y_em,C1,C2,C3,iter] = HDMR_refit(B1,B2,B3,T1,T2,T3,C1,C2, ...
    C3,y_res,Y_em,R,n1,n2,n,m1,m2,m3,j1,j2,j3,ind,vartol,maxiter,maxorder)
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
% Set to zero insignificant components in coefficient matrix + Y_em and   %
% redo backfitting                                                        %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Y_em,C1,C2,C3,iter] = HDMR_refit(B1,B2,B3,T1,T2,T3,C1,C2,C3,...      %
%       y_res,Y_em,R,n1,n2,n,m1,m2,m3,j1,j2,j3,ind,vartol,maxiter,...     %
%       maxorder)                                                         %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Initialize important variables used for convergence analysis backfitting
varb_old = nan(1,n); varb_max = 1; iter = 0; i2 = []; i3 = [];

i1 = find(ind(j1) == 0);                    % indx insignfcnt 1st ord terms
C1(1:m1,i1) = 0; Y_em(1:R,j1(i1)) = 0;      % set id columns to zero
varb_old(1,j1) = sum(C1.^2,1);              % compute old sum coefficients

if maxorder > 1
    i2 = find(ind(j2) == 0);                % indx insignfcnt 2nd ord terms
    C2(1:m2,i2) = 0; Y_em(1:R,j2(i2)) = 0;  % set id columns to zero
    varb_old(1,j2) = sum(C2.^2,1);          % compute old sum coefficients
end

if maxorder > 2
    i3 = find(ind(j3) == 0);                % indx insignfcnt 3rd ord terms
    C3(1:m3,i3) = 0; Y_em(1:R,j3(i3)) = 0;  % set id columns to zero
    varb_old(1,j3) = sum(C3.^2,1);          % compute old sum coefficients
end

% Now collect all coefficients with insignificant sensitivity
ii = [ i1 ; j2(i2) ; j3(i3) ];
% Now determine index of all significant terms
id = 1:n; id(ii) = [];

% Backfitting of all orders combined
while ( varb_max > vartol ) && ( iter < maxiter )
    varb_new = zeros(1,n);
    for i = 1:numel(id)
        % Get column and index of all other non-zero columns
        j = id(i); ii = id; ii(i) = [];
        % Now remove all non-zero columns (nz) from Y_em except column j
        y_nzl(1:R,1) = y_res(1:R,1) - sum(Y_em(1:R,ii),2);
        if ismember(j,j1)
            C1(1:m1,j) = T1(1:m1,1:R,j) * y_nzl(1:R,1);
            Y_em(1:R,j) = B1(1:R,1:m1,j) * C1(1:m1,j);
            varb_new(j) = sum(C1(1:m1,j).^2);
        elseif ismember(j,j2)
            z = j - n1;
            C2(1:m2,z) = T2(1:m2,1:R,z) * y_nzl(1:R,1);
            Y_em(1:R,j) = B2(1:R,1:m2,z) * C2(1:m2,z);
            varb_new(j) = sum(C2(1:m2,z).^2);
        else
            z = j - n1 - n2;
            C3(1:m3,z) = T3(1:m3,1:R,z) * y_nzl(1:R,1);
            Y_em(1:R,j) = B3(1:R,1:m3,z) * C3(1:m3,z);
            varb_new(j) = sum(C3(1:m3,z).^2);
        end
    end
    varb_max = max(abs(varb_new - varb_old));
    varb_old = varb_new;
    iter = iter + 1;
end

end
