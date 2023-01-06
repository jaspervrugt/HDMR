%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  HHH   HHH DDDDDDDD  MMM    MMM RRRRRRR    EEEEEEEE XXX   XXX TTTTTTTTT %
%  HHH   HHH DDDDDDDDD MMM    MMM RRRRRRRR   EEEEEEEE XXX   XXX TTTTTTTTT %
%  HHH   HHH DDD   DDD MMM    MMM RRR   RRR  EEE       XXX XXX  TT TTT TT %
%  HHH   HHH DDD   DDD MMMM  MMMM RRR   RRR  EEE       XXX XXX  T  TTT  T %
%  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRRR   EEEEEE     XXXXX      TTT    %
%  HHHHHHHHH DDD   DDD MMMMMMMMMM RRRRRRR    EEEEEE     XXXXX      TTT    %
%  HHH   HHH DDD   DDD MMM    MMM RRRRRRR    EEE       XXX XXX     TTT    %
%  HHH   HHH DDD   DDD MMM    MMM RRR  RRR   EEE       XXX XXX     TTT    %
%  HHH   HHH DDDDDDDDD MMM    MMM RRR   RRR  EEEEEEEE XXX   XXX    TTT    %
%  HHH   HHH DDDDDDDD  MMM    MMM RRR   RRR  EEEEEEEE XXX   XXX    TTT    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Extended Bases High-Dimensional Model Representation (HDMR) uses three  %
% different type of orthonormal polynomials (Classical, Legendre, and     %
% Hermite) for variance-based global sensitivity analysis (GSA) with      %
% correlated and uncorrelated inputs. This function uses N x d matrix of  %
% N different vector d-vectors of model inputs (factors/parameters) and a %
% N x 1 vector of corresponding model outputs and returns to the user     %
% each factor's first, second, and third order sensitivity coefficient    %
% (separated in total, structural and correlative contributions), an      %
% estimate of their 95% confidence intervals (from bootstrap method)      %
% and the coefficients of the extended bases functions that govern output,%
% Y emulator (determined by an F-test of the error residuals of  the HDMR %
% model with/without a given first, second and/or third order components).%
% These coefficients define an emulator that can be used to predict the   %
% output, Y, of the original model for any d-vector of model inputs. For  %
% uncorrelated model inputs (columns of X are independent), the HDMR      %
% sensitivity indices reduce to a single index (= structural contribution)%
% consistent with their values derived from commonly used variance-based  %
% GSA methods.                                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   SYNOPSIS                                                              %
%    [S,Ss,Fx,Em,XY,RT] = HDMR_EXT(X,Y);                                  %
%    [S,Ss,Fx,Em,XY,RT] = HDMR_EXT(X,Y,options);                          %
%                                                                         %
%   INPUT ARGUMENTS                                                       %
%    X           Nxd matrix: N vectors of d parameters                    %
%    Y           Nx1 vector: single model output for each row of X        %
%    options     (optional) structure: Fields specify HDMR variables      %
%     .graphics  integer [0,1]: graphical output?             (def: 1)    %
%     .basis     integer [1-3]: type of the orthonormal basis (def: 1)    %
%                         1: orthonormal polynomial                       % 
%                          2: legendre polynomial                         % 
%                           3: hermite polynomial                         % 
%     .maxorder  integer [1-3]: max order emulator            (def: 3)    %
%     .m         integer [1-12]: polynomial degree            (def: 3)    %
%     .K         integer [1-500] # bootstrap iterations       (def: 100)  %
%     .R         integer [100-N/2] # bootstrap samples        (def: N/2)  %
%     .alfa      real [0.5-1]: confidence interval F-test     (def: 0.99) %
%     .method    integer [1,2]: 1=forw. sel.; 2=backw. elim.  (def: 1)    %
%                                                                         %
%    def_options = struct('graphics',1,'basis','1','maxorder',3,...       %
%                    'm',3,'K',1,'R',N/2,'alfa',0.99,'method',1);         %
%                                                                         %
%   OUTPUT ARGUMENTS                                                      %
%    S          Cell array: structural, correlative and total sensitivity %
%               each component function -> 1st, 2nd and 3rd order effects %
%               (Table 2 on screen and printed in "HDMR_EXT_results.txt)  %
%    Ss         Cell array: As "S" but lists only significant component   %
%               functions determined via model selection using a F-test   %
%               ( = Table 3 on screen and printed in "HDMR_results.txt" ) %
%    Fx         Cell array: Tabulates emulator properties and performance %
%               on training data set for each of the K bootstrap trials   %
%               ( = Table 1 on screen and printed in "HDMR_results.txt" ) %
%    Em         Structure array: Fields with input and output K emulators %
%     .EB       Nxk matrix: extended bases evaluated at X                 %
%     .C        kx1 matrix: all order coefficients (via D-Morph)          %
%     .m        scalar: degree of orthonormal polynomials                 %
%     .Y_e      RxK matrix: emulator predictions of K bootstraps          %
%     .RMSE     1xK vector: RMSE of emulator residuals of K bootstraps    %
%     .m1       1x1 scalar: number of 1st order terms ( = m )             %
%     .m2       1x1 scalar: number of 2nd order terms ( = 2*m+m^2 )       %
%     .m3       1x1 scalar: number of 3rd order terms ( = 3*m+3*m^2+m^3 ) %
%     .f0       1xK vector: mean Y of each of the K bootstrap trials      %
%     .n        scalar: ( = n1+n2+n3 ) total # of component functions     %
%     .n1       scalar: number of 1st order component functions           %
%     .n2       scalar: number of 2nd order component functions           %
%     .n3       scalar: number of 3rd order component functions           %
%     .k        scalar: total number of terms                             %
%     .k1       scalar: total number of 1st order terms                   %
%     .k2       scalar: total number of 2nd order terms                   %
%     .k3       scalar: total number of 3rd order terms                   %
%     .maxorder scalar: maximum order of HDMR emulator                    %
%     .nterms   1xK vector: # significant terms of each of K emulators    %
%     .p0       1xK vector: # parameters of each of the K emulators       %
%     .RT       1xK vector: CPU time (in sec) to construct each emulator%
%     .select   nxK matrix: significant/insignifant terms each K trials   %
%               if Em.select(i,1) = 0 -> term i insignificant in trial 1  %
%               if Em.select(i,4) = 1 -> term i significant in trial 4    %
%    XY         Structure: X/Y samples and bootstrap information          %
%     .R        scalar: number of random samples each bootstrap trial     %
%     .X_n      Nxd matrix: normalized parameter vectors                  %
%               X_n(i,1:d) = (X(i,1:d)-X_min)./(X_max-X_min); i = 1,..,N  %
%     .minX     1xd vector: min values each X column ( = input )          %
%     .maxX     1xd vector: max values each X column ( = input )          %
%     .Y        Nx1 vector: Y values supplied by user                     %
%     .id       RxK matrix: index R samples of X for each bootstrap trial %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   MAIN REFERENCE                                                        %
%    Genyuan Li, and H. Rabitz, "General Formulation of HDMR component    %
%       functions with independent and correlated variables", Journal of  %
%       Mathematical Chemistry, Vol. 50 Issue 1, pp. 99 - 130, 2012       %
%                                                                         %
%   MATLAB CODE                                                           %
%    Written by Jasper A. Vrugt & Abdullah Sahin                          %
%                                                                         %
%    VERSION: 1.0   FINAL PROGRAM (January 2018)                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First setup HDMR_EXT toolbox (add to path)
addpath(pwd,[pwd,'\miscellaneous']); 
% Then go to an example directory, say example_1
cd example_1
% Now one can execute this example by typing: example_1_uncorrelated