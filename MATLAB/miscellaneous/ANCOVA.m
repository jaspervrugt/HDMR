function [S,S_a,S_b,V_em] = ANCOVA(y,Y_em,V_y,R,n)
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
% ANALYSIS OF COVARIANCE: Returns sensitivity indices each emulator term  %
%                                                                         %
%  SYNOPSIS                                                               %
%   [S,S_a,S_b,V_em] = ANCOVA(y,Y_em,V_Y,R,n)                             %
%  where                                                                  %
%   y           Rx1 vector: single model output for each row of X         %
%   Y_em        Rxn matrix with HDMR model terms from backfitting         %
%   V_y         scalar with variance of the model output                  %
%   R           Number of samples of X                                    %
%   n           Number of terms of HDMR model                             %
%  output                                                                 %
%   S           nx1 vector of total sensitivity each emulator term        %
%   S_a         nx1 vector of structural sensitivity each emulator term   % 
%   S_b         nx1 vector of correlative sensitivity each emulator term  % 
%   V_em        nx1 vector of variance each emulator term                 %
%                                                                         %
%  REFERENCE                                                              %
%   Genyuan Li et al., Journal of Physical Chemistry A., V. 114 (19),     %
%       pp. 6022-6032, 2010.                                              %
%                                                                         %
%  © Written by Jasper A. Vrugt, Sept. 9, 2017                            %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% NOTE: f0 excluded from Y_bf and Y_ind as covariance unaffected by constant

% Now compute the sum of all Y_bf terms
Y0 = sum(Y_em(1:R,1:n),2);
% Initialize each variable
[S,S_a,S_b,V_em] = deal(nan(n,1));

% Analysis of covariance -> extension of analysis of variance
for j = 1:n
    C = cov(Y_em(1:R,j),y);                 % Covariance matrix of jth term of Y_bf and actual Y 
    S(j,1) = C(1,2)/V_y;                    % Total sensitivity of jth term         ( = Eq. 19 of Li et al )
    C = cov(Y_em(1:R,j),Y0 - Y_em(1:R,j));  % Covariance matrix of jth term with emulator Y without jth term
    S_a(j,1) = C(1,1)/V_y;                  % Structural contribution of jth term   ( = Eq. 20 of Li et al ) 
    S_b(j,1) = C(1,2)/V_y;                  % Correlative contribution of jth term  ( = Eq. 21 of Li et al )
    V_em(j,1) = sum(C(1,1:2));              % Variance in Y of jth term             ( = S_a * V_Y + S_b * V_Y )
                                            %                                       ( = (S_a + S_b) * V_Y )
                                            %                                       ( = S * V_Y = C(1,2) L37 )
end

end
