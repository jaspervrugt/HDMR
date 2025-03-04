function [ind,nterms,p0] = HDMR_F_test(Y,f0,Y_bf,R,alfa,m1,m2,m3,n1,n2, ...
    n3,n,method)
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
% Model selection using the F test                                        %
%                                                                         %
%  SYNOPSIS                                                               %
%   [ind,nterms,p0] = HDMR_F_test(Y,f0,Y_bf,R,alfa,m1,m2,m3,n1,n2,n3, ... %
%       n,method)                                                         %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Initialize ind with zeros (all terms insignificant)
ind = zeros(n,1);

% Determine the significant components of the HDMR model via the F-test
switch method
    case 1 % --> forward selection of terms (start with Y = f0)
        Y_res0 = Y - f0; SSR0 = sum ( Y_res0.^2 ); p0 = 0;
        for i = 1:n
            % model with ith term included
            Y_res1 = Y_res0 - Y_bf(:,i);
            % Number of parameters of proposed model (order dependent)
            if i <= n1           
                p1 = p0 + m1;        % 1st order
            elseif i > n1 && i <= n1 + n2
                p1 = p0 + m2;       % 2nd order
            else
                p1 = p0 + m3;       % 3rd order
            end
            % Calculate SSR of Y1
            SSR1 = sum( ( Y_res1 ).^2 ); 
            % Now calculate the F_stat (F_stat > 0 -> SSR1 < SSR0 )
            F_stat = ( (SSR0 - SSR1)/(p1 - p0) ) / ( SSR1/(R - p1) );
            % Now calculate critical F value
            F_crit = finv(1 - alfa,p1 - p0,R - p1);
            % Now determine whether to accept ith component into model
            if F_stat > F_crit
                % ith term is significant and should be included in model
                ind(i) = 1; Y_res0 = Y_res1; SSR0 = SSR1; p0 = p1;
            end
        end
    case 2 % --> backward elimination of terms (start with Y = f0 + sum(all_terms))
        % NOTE: ONLY POSSIBLE IF R - p1 > 0 OTHERWISE MUST DO FORWARD SELECTION
        % Determine output of full model (all terms included)
        Y_res0 = Y - f0 - sum(Y_bf,2); 
        % Calculate the SSR of the full model (all terms included)
        SSR0 = sum ( Y_res0.^2 );
        % Determine number of parameters of full model
        p0 = n1 * m1 + n2 * m2 + n3 * m3;
        for i = 1:n
            % previous model with ith term excluded
            Y_res1 = Y_res0 + Y_bf(:,i);
            % Number of parameters of proposed model (order dependent)
            if i <= n1
                p1 = p0 - m1;           % 1st order
            elseif i > n1 && i <= n1 + n2
                p1 = p0 - m2;           % 2nd order
            else
                p1 = p0 - m3;           % 3rd order        
            end
            % Calculate SSR of Y1
            SSR1 = sum( Y_res1.^2 );
            % Now calculate the F_stat (F_stat > 0 if SS1 > SSR0 as p1 < p0)
            F_stat = ( (SSR0 - SSR1)/(p1 - p0) ) / ( SSR1/(R - p1) );
            % Now calculate critical F value
            F_crit = finv(1 - alfa,p0 - p1,R - p1); % Note 2nd term turned around
            % Now determine whether ith component is significant
            if F_stat > F_crit
                % ith term is significant and should retain in model!
                ind(i) = 1;
            else
                % ith term is insignificant and will be removed from model
                p0 = p1; SSR0 = SSR1; Y_res0 = Y_res1; 
            end
        end
end

% Now return number of terms of final HDMR model
nterms = sum(ind);