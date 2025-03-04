function [SA,SA_sig,Fx,Em,Xy,RT] = HDMR(X,y,options)
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
%                                                                         %
%  High-Dimensional Model Representation (HDMR) using B-spline functions  %
%  for variance-based global sensitivity analysis (GSA) with correlated   %
%  and uncorrelated inputs. This function uses as input a N x d matrix    %
%  of N different d-vectors of model inputs (factors/parameters) and a    %
%  N x 1 vector of corresponding model outputs and returns to the user    %
%  each factor's first, second, and third order sensitivity coefficient   %
%  (separated in total, structural and correlative contributions), an     %
%  estimate of their 95% confidence intervals (from bootstrap method)     %
%  and the coefficients of the significant B-spline basis functions that  %
%  govern output, y (determined by an F-test of the error residuals of    %
%  the HDMR model (emulator) with/without a given first, second and/or    %
%  third order B-spline). These coefficients define an emulator that can  %
%  be used to predict the output, y, of the original (CPU-intensive?)     %
%  model for any d-vector of model inputs. For uncorrelated model inputs  %
%  (columns of X are independent), the HDMR sensitivity indices reduce    %
%  to a single index (= structural contribution), consistent with their   %
%  values derived from commonly used variance-based GSA methods.          %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  SYNOPSIS                                                               %
%   [S,Ss,Fx,Em,Xy,RT] = HDMR(X,y);                                       %
%   [S,Ss,Fx,Em,Xy,RT] = HDMR(X,y,options);                               %
%  where                                                                  %
%   X           [input] Nxd matrix: N vectors of d parameters             %
%   y           [input] Nx1 vector: single model output each row X        %
%   options     [input] (optional) structure: HDMR variables              %
%    .graphics   integer [0,1]: graphical output?             (def: 1)    %
%    .maxorder   integer [1-3]: max order emulator            (def: 3)    %
%    .maxiter    integer [1-1000]: max iterations backfitting (def: 100)  %
%    .bf1        integer [0,1]: 1st order, ind./backfitting   (def: 1)    %
%    .bf2        integer [0,1]: 2nd order, ind./backfitting   (def: 1)    %
%    .bf3        integer [0,1]: 3rd order, ind./backfitting   (def: 1)    %
%    .method     integer [1,2]: 1=forw. sel.; 2=backw. elim.  (def: 1)    %
%    .m          integer [2-10]: # B-spline intervals         (def: 2)    %
%    .K          integer [1-500] # bootstrap iterations       (def: 100)  %
%    .R          integer [100-N/2] # bootstrap samples        (def: N/2)  %
%    .alfa       real [0.-0.1]: significance level F-test     (def: 0.01) %
%    .lambda     real [0-inf]: regularization coefficient     (def: 0.10) %
%    .vartol     real (0-1]: tolerance backfitting            (def: 1e-5) %
%   def_options = struct('graphics',1,'maxorder',3,'maxiter','100',...    %
%                   'bf1',1,'bf2',1,'bf3',1,'method',1,'m',2,'K',100,...  %
%                   'R',N/2,'alfa',0.01,'lambda',0.1,'vartol',1e-5);      %
%   S           [outpt] Cell array: structural, correlative & total       %
%                       sensitivity each component function -> 1st, 2nd   %
%                       and 3rd order effects                             %
%                       (= Table 2 screen and in "HDMR_results.txt")      %
%   Ss          [outpt] Cell array: As "S" but lists only significant     %
%                       component functions determined via model          %
%                       selection using a F-test                          %
%                       (= Table 3 screen and in "HDMR_results.txt")      %
%   Fx          [outpt] Cell array: Tabulates emulator properties and     %
%                       performance on training data set for each of the  %
%                       K bootstrap trials                                %
%                       (= Table 1 screen and in "HDMR_results.txt" )     %
%   Em          [outpt] Structure array: Fields input/output K emulators  %
%    .B1         Nx(m+3)xn1 matrix: B-spline function evaluated at X      %
%    .B2         Nx(m+3)^2xn2 matrix: 2nd order B-spline evaluated at X   %
%    .B3         Nx(m+3)^3xn3 matrix: 3rd order B-spline evaluated at X   %
%    .C1         dxn1xK matrix: 1st order coefficients (backfitting)      %
%    .C2         d2xn2xK matrix: 2nd order coefficients (backfiting)      %
%    .C3         d^3xn3xK matrix: 3rd order coefficients (backfitting)    %
%    .m          scalar: # spline intervals ( =  + 3)                     %
%    .Y_e        RxK matrix: Emulator predictions K training data sets    %
%    .RMSE       1xK vector: RMSE of emulator residuals K training sets   %
%    .c1         n1x1 vector: indices 1st order terms ( = 1:d )           %
%    .c2         n2x2 matrix: indices 2nd order combinations              %
%    .c3         n3x3 matrix: indices 3rd order combinations              %
%    .f0         1xK vector: mean y of each of K bootstrap trials         %
%    .n          scalar: ( = n1+n2+n3 ) max. number component functions   %
%    .n1         scalar: ( = d) with total number of 1st order terms      %
%    .n2         scalar: total number of 2nd order component functions    %
%    .n3         scalar: total number of 3rd order component functions    %
%    .maxorder   scalar: maximum order of HDMR emulator                   %
%    .nterms     1xK vector: # significant terms of each of K emulators   %
%    .p0         1xK vector: # parameters of each of the K emulators      %
%    .select     nxK matrix: significant/insignifant terms each K trials  %
%                  if Em.select(i,1) = 0 -> term i insignificant trial 1  %
%                  if Em.select(i,4) = 1 -> term i significant trial 4    %
%   Xy          [outpt] Structure: X/y samples and bootstrap information  %
%    .R          scalar: number of random samples each bootstrap trial    %
%    .X_n        Nxd matrix: normalized parameter vectors                 %
%                X_n(i,1:d) = (X(i,1:d)-X_min)./(X_max-X_min); i = 1,..,N %
%    .X_min      1xd vector: min values each X column ( = input )         %
%    .X_max      1xd vector: max values each X column ( = input )         %
%    .y          Nx1 vector: Y values supplied by user                    %
%    .id         RxK matrix: index R samples of X each bootstrap trial    %
%   RT          [outpt] 1xK vector: CPU time (sec) to construct emulator  %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  BUILT-IN CASE STUDIES                                                  %
%   Example 1   Multivariate normal benchmark study                       %
%   Example 2   Multivariate normal benchmark study with Σ ≠ identity     %
%   Example 3   Multivariate normal benchmark study with Σ ≠ identity     %
%   Example 4   Ishigami function                                         %
%   Example 5   Ishigami function with Σ ≠ identity                       %
%   Example 6   Function with Σ ≠ identity                                %
%   Example 7   Sobol g function                                          %
%   Example 8   Multivariate normal with Σ ≠ identity                     %
%   Example 9   Example 1 of Chastaing et al. (2012)                      %
%   Example 10  Example 2 of Chastaing et al. (2012)                      %
%   Example 21  Soil temperature modeling                                 %
%   Example 22  Rainfall runoff modeling: hmodel                          %
%   Example 23  Rainfall runoff modeling: SAC-SMA model                   %
%   Example 24  Dynamic foodweb: One-predator-one-prey model              %
%   Example 25  Dynamic foodweb: Two-predators-two-preys model            %
%   Example 26  Dynamic foodweb: Two-predators-two-preys model real data  %
%   Example 27  Simple multivariate function                              %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  MAIN REFERENCE                                                         %
%   Gao, Y., C. Guilloteau, E. Foufoula-Georgiou, C. Xu, X. Sun, and      %
%       J.A. Vrugt (2024), Soil moisture-cloud-precipitation feedback in  %
%       the lower atmosphere from functional decomposition of satellite   %
%       observations, Geophysical Research Letters, 51, e2024GL110347,    %
%       https://doi.org/10.1029/2024GL110347                              %
%   Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic sensitivity   %
%       analysis with dependent variables: Covariance-based decomposition %
%       of hydrologic models, Water Resources Research, 59 (4),           %
%       e2022WR0328346, https://doi.org/10.1029/2022WR032834              %
%   Li, G. H. Rabitz, P.E. Yelvington, O.O. Oluwole, F. Bacon,            %
%       C.E. Kolb, and J. Schoendorf (2010), Global sensitivity analysis  %
%       for systems with independent and/or correlated inputs, Journal of %
%       Physical Chemistry A, Vol. 114 (19), pp. 6022 - 6032              %
%                                                                         %
%  MATLAB CODE                                                            %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
%  Version 1.0    January 2018                                            %
%  Version 2.0    July 2024                                               %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 2
    error(['HDMR ERROR: Insufficient number of input arguments. ' ...
        'Use "help install_HDMR"']);
else
    close all hidden;
end

% INITIAL PART: Check arguments & define HDMR variables
y_in = y(:); ct = 0;                                    % Make sure that Y is a vertical vector and ct is zero
[N,d,graphics,maxorder,maxiter,bf1,bf2,bf3,m,K,R,...    % Check the values/content of structure options    
    method,alfa,vartol,lambda,refit] = ...
    HDMR_setup(X,y_in,options);  
[Xy,Em,SA,RT,Y_em,T2,T3,m1,m2,m3,j1,j2,j3,it2,...       % Initialize main variables used by HDMR
    it3,itr] = HDMR_initialize(X,y_in,N,d,K,R,m,...
    maxorder);         

for k = 1:K     % DYNAMIC PART: Bootstrap for confidence intervals

    start_time = cputime;                               % Start timer   
    y = y_in(Xy.id(1:R,k),1);                           % Extract the "right" y values
    SA.V_y(k) = var(y);                                 % Compute the variance of the model output
    Em.f0(k) = sum(y)/R;                                % Mean of the output
    Y_res = y - Em.f0(k);                               % Compute residual
    [Y_res,Y_em(1:R,j1),Em.C1(:,:,k),T1,it1] = ...      % 1st order component functions: ind/backfitting
        HDMR_1st(Em.B1(Xy.id(1:R,k),1:m1,1:Em.n1),...
        Y_res,R,Em.n1,m1,vartol,lambda,maxiter,bf1);    
    if ( maxorder > 1 )
        [Y_res,Y_em(1:R,j2),Em.C2(:,:,k),T2,it2] = ...  % 2nd order component functions: ind/backfitting
        HDMR_2nd(Em.B2(Xy.id(1:R,k),:,:),Y_res,...
            R,Em.n2,m2,vartol,lambda,maxiter,bf2);      
    end
    if ( maxorder == 3 )
        [Y_em(1:R,j3),Em.C3(:,:,k),T3,it3] = ...        % 3rd order component functions: ind/backfitting
            HDMR_3rd(Em.B3(Xy.id(1:R,k),:,:),Y_res,...
            R,Em.n3,m3,vartol,lambda,maxiter,bf3);      
    end
    [Em.select(1:Em.n,k),Em.nterms(k),Em.p0(k)] = ...   % Identify significant and insignificant terms
        HDMR_F_test(y,Em.f0(k),Y_em,R,alfa,m1,m2, ...
        m3,Em.n1,Em.n2,Em.n3,Em.n,method);    
    if ( refit == 1 )
        [Y_em,Em.C1(:,:,k),Em.C2(:,:,k), ...            % Do simultaneous backfitting of all terms
            Em.C3(:,:,k),itr] = HDMR_refit( ...
            Em.B1(Xy.id(1:R,k),1:m1,1:Em.n1), ...
            Em.B2(Xy.id(1:R,k),:,:), ...
            Em.B3(Xy.id(1:R,k),:,:),T1,T2,T3, ...
            Em.C1(:,:,k),Em.C2(:,:,k), ...
            Em.C3(:,:,k),y - Em.f0(k),Y_em,R, ...
            Em.n1,Em.n2,Em.n,m1,m2,m3,j1,j2,j3, ...
            Em.select(1:Em.n,k),vartol,maxiter, ...
            maxorder);                   
    end
    Em.Y_e(:,k) = Em.f0(k) + sum(Y_em,2);               % Store the emulated output
    Em.RMSE(k) = sqrt(sum((y - Em.Y_e(1:R,k)).^2)/R);   % RMSE of emulator
    RT(k) = cputime - start_time;                       % CPU time kth emulator
    Em.iter(1:4,k) = [ it1 ; it2 ; it3 ; itr ];         % # iterations of 1st, 2nd, 3rd, and refit
    [SA.S(1:Em.n,k),SA.Sa(1:Em.n,k), ...                % Compute sensitivity indices
        SA.Sb(1:Em.n,k),SA.V_em(1:Em.n,k)] = ...
        ANCOVA(y,Y_em,SA.V_y(k),R,Em.n);                
    fprintf(1, repmat('\b',1,ct));                      % Print progress
    ct = fprintf(['HDMR CALCULATING: FINISHED TRIAL ' ...
        '%d OF %d AND THUS %3.2f%% DONE'],k,K,100*(k/K));
    
end         % DYNAMIC PART: End of bootstrap

% FINAL PART: Postprocess in tables & figures
[SA,SA_sig,Fx] = HDMR_end(SA,Em.nterms,Em.p0, ...       % Now call function HDMR_end
    Em.RMSE,Em.select,K,Em.c2,Em.c3,Em.n1, ...
    Em.n2,Em.n3,Em.n);               
if graphics == 1
    HDMR_plot(SA_sig,Fx,y_in,Em.Y_e,Em.select,...       % If graphical output is desired: plotting function
        Em.p0,Em.iter,Xy.id,R,K,refit,maxorder);        
end

end
