function [Xy,Em,SA,RT,Y_em,T2,T3,m1,m2,m3,j1,j2,j3,it2,it3,itr] = ...
    HDMR_initialize(X,y,N,d,K,R,m,maxorder)
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
% Initialize main variables used by HDMR                                  %
%                                                                         %
%  SYNOPSIS                                                               %
%   [Xy,Em,SA,RT,Y_em,T2,T3,m1,m2,m3,j1,j2,j3,it2,it3,itr] = ...          %
%       HDMR_initialize(X,y,N,d,K,R,m,maxorder)                           %
%  where                                                                  %
%    X           Nxd matrix: N vectors of d parameters                    %
%    y           Nx1 vector: single model output for each row of X        %
%                                                                         %
%  REFERENCE                                                              %
%   Genyuan Li, H. Rabitz, P.E. Yelvington, O.O. Oluwole, F. Bacon,       %
%       C.E. Kolb, and J. Schoendorf, "Global Sensitivity Analysis for    %
%       Systems with Independent and/or Correlated Inputs", Journal of    %
%       Physical Chemistry A, Vol. 114 (19), pp. 6022 - 6032, 2010        %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Random seed (legacy: randn('state',sum(100*clock)); )
rng(1+round(100*rand),'twister');

% STRUCTURE XY: Define content
if K == 1
    Xy = struct('X_n',nan(N,d),'minX',min(X,[],1),'maxX', ...
        max(X,[],1),'y',y,'R',R,'id',(1:R)');
else
    % Now setup the boostrap matrix with selection matrix, id, for samples
    [~,id] = sort(rand(N,K));
    % Store in structure XY
    Xy = struct('X_n',nan(N,d),'minX',min(X,[],1),'maxX', ...
        max(X,[],1),'y',y,'R',R,'id',id(1:R,1:K));
end
% Compute normalized X-values
Xy.X_n = (X(1:N,1:d) - repmat(Xy.minX(1:d),N,1)) ./ ...
    repmat(Xy.maxX(1:d) - Xy.minX(1:d),N,1);

% STRUCTURE Em: Important variables
[n2,n3] = deal(0); [c2,c3] = deal([]);
% determine all combinations of parameters (coefficients) for each order
c1 = 1:d; n1 = d; % = size(c1,1);
% now return based on maxorder used
if maxorder > 1
    c2 = nchoosek(1:d,2); n2 = size(c2,1);
end
if maxorder == 3
    c3 = nchoosek(1:d,3); n3 = size(c3,1);
end
% calulate total number of coefficients
n = n1 + n2 + n3;

% Initialize m1, m2 and m3 - number of coefficients first, second, third order
m1 = m + 3; m2 = m1^2; m3 = m1^3;

% STRUCTURE Em: Initialization
switch maxorder
    case 1
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE', ...
            nan(1,K),'m',m,'Y_e',nan(R,K),'f0',nan(1,K),...
            'c1',c1,'n1',d,'c2',c2,'n2',n2,'c3',c3,'n3',n3,'n',n, ...
            'maxorder',maxorder,'select',nan(n,K),...
            'C1',nan(m1,n1,K),'C2',nan(1,1,K),'C3',nan(1,1,K),'B1', ...
            zeros(N,m1,n1),'B2',...
            nan(N,1),'B3',nan(N,1),'iter',nan(4,K));
    case 2
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),...
            'c1',c1,'n1',n1,'c2',c2,'n2',n2,'c3',c3,'n3',n3,'n',n, ...
            'maxorder',maxorder,'select',nan(n,K),...
            'C1',nan(m1,n1,K),'C2',nan(m2,n2,K),'C3',nan(1,1,K), ...
            'B1',zeros(N,m1,n1),'B2',...
            zeros(N,m2,n2),'B3',nan(N,1),'iter',nan(4,K));
    case 3
        Em = struct('nterms',nan(1,K),'p0',nan(1,K),'RMSE',nan(1,K), ...
            'm',m,'Y_e',nan(R,K),'f0',nan(1,K),...
            'c1',c1,'n1',n1,'c2',c2,'n2',n2,'c3',c3,'n3',n3,'n',n, ...
            'maxorder',maxorder,'select',nan(n,K),...
            'C1',nan(m1,n1,K),'C2',nan(m2,n2,K),'C3',nan(m3,n3,K), ...
            'B1',zeros(N,m1,n1),'B2',...
            zeros(N,m2,n2),'B3',zeros(N,m3,n3),'iter',nan(4,K));
end

% Now compute B-spline values for all N samples of X_n
Em.B1 = B_spline(Xy.X_n,N,d,m); 

% Now compute B values for second order
if maxorder > 1
    beta = permn(1:m1,2);
    for k = 1:n2
        for j = 1:m2
            Em.B2(1:N,j,k) = Em.B1(1:N,beta(j,1),Em.c2(k,1)) ...
                .* Em.B1(1:N,beta(j,2),Em.c2(k,2));
        end
    end
end

% Compute B values for third order
if maxorder == 3
    % Now compute B values for third order
    beta = permn(1:m1,3);
    for k = 1:n3
        for j = 1:m3
            Em.B3(:,j,k) = Em.B1(:,beta(j,1),Em.c3(k,1)) ...
                .* Em.B1(:,beta(j,2),Em.c3(k,2)) .* ...
                Em.B1(:,beta(j,3),Em.c3(k,3));
        end
    end
end

% STRUCTURE SA: Sensitivity analysis and analysis of variance decomposition
SA = struct('S',nan(Em.n,K),'Sa',nan(Em.n,K),'Sb',nan(Em.n,K), ...
    'ST',nan(d,1),'V_em',nan(Em.n,K),'V_Y',nan(1,K));
% Return runtime
RT = zeros(1,K);
% Initialize various terms
Y_em = nan(R,Em.n); [T2,T3] = deal([]);
% Initialize index of first, second and third order terms (columns of Y_bf)
j1 = (1:n1)'; j2 = (n1+1:n1+n2)'; j3 = (n1+n2+1:n)';
% Now initialize number of iterations first, second and third order
[it2,it3,itr] = deal(0);

% Print to screen
fprintf('\n');
fprintf('  ---------------------------------------------------------------------  \n');
fprintf('           HHH     HHH   DDDDDDDDD    MMM     MMM   RRRRRRRR             \n');
fprintf('           HHH     HHH   DDDDDDDDDD   MMM     MMM   RRRRRRRRR            \n');
fprintf('           HHH     HHH   DDD    DDD   MMMM   MMMM   RRR    RRR           \n');
fprintf('           HHH     HHH   DDD    DDD   MMMMM MMMMM   RRR    RRR           \n');
fprintf('           HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRRR            \n');
fprintf('           HHHHHHHHHHH   DDD    DDD   MMMMMMMMMMM   RRRRRRRR             \n');
fprintf('           HHH     HHH   DDD    DDD   MMM     MMM   RRRRRRRRR            \n');
fprintf('           HHH     HHH   DDD    DDD   MMM     MMM   RRR   RRR            \n');
fprintf('           HHH     HHH   DDDDDDDDDD   MMM     MMM   RRR    RRR           \n');
fprintf('           HHH     HHH   DDDDDDDDD    MMM     MMM   RRR     RRR          \n');
fprintf('  ---------------------------------------------------------------------  \n');
fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
fprintf('\n'); 

% % % Print to screen
% % fprintf('\n');
% % fprintf(['  ---------------------------------------------------------' ...
% %     '--------------------------------------------------  \n']);
% % fprintf('      H     H  DDDDD    M     M  RRRRRR     PPPP        A     CCCCCCC  K      KK     A     GGGGGG  EEEEEE \n');
% % fprintf('      H     H  DDDDDD   M     M  RRRRRRR    PPPPP      AAA    CCCCCC   K     KK     AAA    G    G  EEEEEE \n');  
% % fprintf('      H     H  DDDDDDD  MM   MM  R     R    PPPPPP    AAAAA   CC       K     KK    AAAAA   G    G  E      \n');    
% % fprintf('      H     H  D     D  MMM MMM  R    R     P    PP  A     A  C        K    KK    A     A  G    G  E      \n');  
% % fprintf('      HHHHHHH  D     D  MMMMMMM  R   R      P    PP  A     A  C        KKKKKK     A     A  GGGGGG  EEE    \n');  
% % fprintf('      HHHHHHH  D     D  M MMM M  RRRRR      PPPPPP   AAAAAAA  C        KKKKKK     AAAAAAA  GGGGGG  EEE    \n');    
% % fprintf('      H     H  D     D  M  M  M  RRRRR      P        AAAAAAA  C        K    KK    AAAAAAA       G  E      \n');  
% % fprintf('      H     H  DDDDDDD  M     M  R   R      P        A     A  CC       K     KK   A     A       G  E      \n');  
% % fprintf('      H     H  DDDDDD   M     M  R    R     P        A     A  CCCCCC   K     KK   A     A      GG  EEEEEE \n');    
% % fprintf('      H     H  DDDDD    M     M  R     R    P        A     A  CCCCCCC  K      KK  A     A  GGGGGG  EEEEEE \n');    
% % fprintf(['  ---------------------------------------------------------' ...
% %     '--------------------------------------------------  \n']);
% % fprintf('  © Jasper A. Vrugt, University of California Irvine \n');
% % fprintf('\n'); 

end
