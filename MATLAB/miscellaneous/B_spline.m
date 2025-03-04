function B = B_spline(X,R,d,m)
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
% Returns coefficients, B, of B-splines. Written to keep in mind C++      %
% translation as vectorization can make things faster in MATLAB           %
%                                                                         %
%  SYNOPSIS                                                               %
%   B = B_spline(X,R,d,m)                                                 %
%  where                                                                  %
%   X           Rxd matrix: R vectors of d parameters                     %
%   R           # samples of X                                            %
%   d           # parameters                                              %
%   m           # B-spline intervals                                      %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Initialize matrix B 
B = zeros(R,m+3,d); % B = nan(R,M+3,d); % --> uncomment Line 17 - 19
% Calculate interval, h1
h1 = 1/m;
% Now loop over each parameter of X and B-spline interval
for i = 1:d
    for j = 1:m+3
        % Change value of k
        k = j - 2;
        for r = 1:R
            % Check interval of X(r,i)
%             if X(r,i) <= ( k - 2 ) * h1 || X(r,i) > ( k + 2 ) * h1
%                 B(r,m,i) = 0;
%             end 
            % Check another interval of X(r,i)        
            if X(r,i) > ( k + 1 ) * h1 && X(r,i) <= ( k + 2 ) * h1
                B(r,j,i) = ( ( k + 2 ) * h1 - X(r,i) ).^3;
            end
            % Check another interval of X(r,i)        
            if X(r,i) > k * h1 && X(r,i) <= ( k + 1 ) * h1
                B(r,j,i) = ( ( k + 2 ) * h1 - X(r,i) ).^3 - ...
                    4 * ( ( k + 1 ) * h1 - X(r,i) ).^3;
            end
            % Check another interval of X(r,i)        
            if X(r,i) > ( k - 1 ) * h1 && X(r,i) <= k * h1
                B(r,j,i) = ( ( k + 2 ) * h1 - X(r,i) ).^3 - ...
                    4 * ( ( k + 1 ) * h1 - X(r,i) ).^3 + ...
                    6 * ( k * h1 - X(r,i) ).^3;
            end
            % Check another interval of X(r,i)        
            if X(r,i) > ( k - 2 ) * h1 && X(r,i) <= ( k - 1 ) * h1
                B(r,j,i) = ( ( k + 2 ) * h1 - X(r,i) ).^3 - ...
                    4 * ( ( k + 1 ) * h1 - X(r,i) ).^3 + ...
                    6 * ( k * h1 - X(r,i) ).^3 - ...
                    4 * ( ( k - 1 ) * h1 - X(r,i) ).^3;
            end
        end
    end
end
% Multiply B with m^3
B = m.^3 * B;