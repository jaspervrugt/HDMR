function varargout = HDMR_model(Em,Xy,X_evl)
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
% Returns y values of HDMR emulator based on coefficients from HDMR       %
%                                                                         %
%  SYNOPSIS                                                               %
%   [y_e,Y_ebf,y_evl,Y_evlbf] = HDMR_model(Em,Xy,X_evl)                   %
%  where                                                                  %
%    Em         [input] Emulator structure from HDMR package              %
%    Xy         [input] Sample structure from HDMR package                %
%    X_evl      [input] OPT: New samples of X not used by HDMR package    %
%    y_e        [outpt] Nx1 vector HDMR emulator last bootstrap trial     %
%    Y_ebf      [outpt] Nxnc matrix HDMR model terms last bootstrap trial %
%    y_evl      [outpt] OPT: Nx1 vector of emulator predictions of X_evl  %
%    Y_evlbf    [outpt] OPT: NxNc matrix HDMR model terms relatd to X_evl %
%                                                                         %
%  © Written by Jasper A. Vrugt                                           %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Print to screen function help
fprintf('\n');
fprintf('%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< %%\n');
fprintf('%% HDMR_model computes the HDMR emulator of last bootstrap trial           %%\n');
fprintf('%%   with original samples of X and Y stored in structure XY               %%\n');
fprintf('%%       and/or a new set of samples stored in X_new                       %%\n');
fprintf('%%                                                                         %%\n');
fprintf('%%  SYNOPSIS                                                               %%\n');
fprintf('%%   [y_e,Y_ebf,y_evl,Y_evlbf] = HDMR_model(Em,Xy,X_evl)                   %%\n');
fprintf('%%  where                                                                  %%\n');
fprintf('%%    Em         [input] Emulator structure from HDMR package              %%\n');
fprintf('%%    Xy         [input] Sample structure from HDMR package                %%\n');
fprintf('%%    X_evl      [input] OPT: New samples of X not used by HDMR package    %%\n');
fprintf('%%    y_e        [outpt] Nx1 vector HDMR emulator last bootstrap trial     %%\n');
fprintf('%%    Y_ebf      [outpt] Nxnc matrix HDMR model terms last bootstrap trial %%\n');
fprintf('%%    y_evl      [outpt] OPT: Nx1 vector of emulator predictions of X_evl  %%\n');
fprintf('%%    Y_evlbf    [outpt] OPT: NxNc matrix HDMR model terms relatd to X_evl %%\n');
fprintf('%%                                                                         %%\n');
fprintf('%%  © Written by Jasper A. Vrugt                                           %%\n');
fprintf('%%    University of California Irvine                                      %%\n');
fprintf('%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>< %%\n');
fprintf('\n');

% Check number of input arguments
if nargin < 2
    error(['HDMR_EMULATE ERROR: Insufficient number of inputs. ' ...
        'Please use: [y_e,Y_ebf] = HDMR_model(Em,Xy)']);
end

% Unpack XY
X = Xy.X_n; y = Xy.y;
% Determine size of X
[N,d] = size(X);
% Determine K, number of bootstrap trials
K = size(Em.C1,3); R = Xy.R;
% extract id of selected samples last bootstrap trial
id = Xy.id(:,K);
% Unpack the emulator coefficients from Fx
if isfield(Em,'C1b')
    [C1,C2,C3] = deal(Em.C1b(:,:,K),Em.C2b(:,:,K),Em.C3b(:,:,K));
else
    [C1,C2,C3] = deal(Em.C1(:,:,K),Em.C2(:,:,K),Em.C3(:,:,K));
end
% Set f0 to be mean of Y of last emulator
f0 = Em.f0(K);
% Now built the last emulator of the bootstrap trials
[y_e,Y_ebf] = HDMR_built(N,d,f0,Em.B1,Em.B2,Em.B3,C1, ...
    C2,C3,Em.n1,Em.n2,Em.n3,Em.n);
% Now create plot with emulator output
HDMR_plot(y_e,y,id,N,R,K);

% Check whether we have new evaluation samples as input
if nargin == 3
    % Evaluate the evaluation samples in X_new
    [N,d] = size(X_evl);
    % Normalize the X values
    X_n(1:N,1:d) = (X_evl(1:N,1:d) - repmat(Xy.minX(1:d),N,1)) ./ ...
        repmat(Xy.maxX(1:d)-Xy.minX(1:d),N,1);
    % Calculate B splines
    B1 = B_spline(X_n,N,d,Em.m);
    % Return function that computes beta for 2nd and 3rd order
    f_beta = @(M,order) dec2base(0:(M+3)^order - 1,M+3,order) - 47;
    % Calculate B2
    beta2 = f_beta(Em.m,2); B2 = zeros(N,(Em.m+3)^2,Em.n2);
    for k = 1:Em.n2
        for j = 1:(Em.m+3)^2
            B2(1:N,j,k) = B1(1:N,beta2(j,1),Em.c2(k,1)) .* ...
                B1(1:N,beta2(j,2),Em.c2(k,2));
        end
    end
    % Calculate B3
    beta3 = f_beta(Em.m,3); B3 = zeros(N,(Em.m+3)^3,Em.n3);
    for k = 1:Em.n3
        for j = 1:(Em.m+3)^3
            B3(1:N,j,k) = B1(1:N,beta3(j,1),Em.c3(k,1)) .* ...
                B1(1:N,beta3(j,2),Em.c3(k,2)) .* ...
                B1(1:N,beta3(j,3),Em.c3(k,3));
        end
    end
    % Set f0 to be mean of Y of last emulator
    f0 = sum(Em.f0)/K;
    % Built the last emulator of the bootstrap trials
    [y2_e,Y2_ebf] = HDMR_built(N,d,f0,B1,B2,B3,C1,C2, ...
        C3,Em.n1,Em.n2,Em.n3,Em.n);
    % Return output arguments
    varargout(1:4) = { y_e Y_ebf y2_e Y2_ebf };
else
    % Return output arguments
    varargout(1:2) = { y_e Y_ebf };
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary functions
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% First: HDMR_built
function [y_e,Y_ebf] = HDMR_built(N,d,f0,B1,B2,B3,coef1, ...
    coef2,coef3,c1,c2,c3,nc)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%HDMR_built: Computes the terms of the emulator
% Written by Jasper A. Vrugt
% University of California Irvine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Initialize the HDMR component functions
Y_ebf = nan(N,nc);

for z = 1:c1        % --> first order terms
    Y_ebf(1:N,z) = B1(:,:,z) * coef1(:,z);
end

for k = 1:c2        % --> second order terms
    for z = d + 1 : d + c2
        Y_ebf(1:N,z) = B2(:,:,z - d) * coef2(:,z - d);
    end
end

for k = 1:c3        % --> third order terms
    for z = d + 1 + c2 : nc
        Y_ebf(1:N,z) = B3(:,:,z - (d + c2)) * coef3(:,z - (d + c2));
    end
end

% Now compute the sum of the Y_ebf values
y_e = f0 + sum(Y_ebf,2);

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% First: HDMR_plot
function HDMR_plot(y_e,y,id,N,R,K)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%HDMR_plot: Computes the terms of the emulator
% Written by Jasper A. Vrugt
% University of California Irvine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Set the default interpreter to latex
set(0,'defaultTextInterpreter','latex');
% Create a plot
figure('units','normalized','outerposition',[0 0 1 1])
% Now generate y values
y_plot = linspace(min(y),max(y));

% Define multiplier to round RMSE values
T = 1000;

% Determine screen size
scr_z = get(0,'ScreenSize'); 
% Multiplier, x and y: axis
x_mult = scr_z(3)/1920; y_mult = scr_z(4)/1080;
% Multiplier, text
t_mult = min(x_mult,y_mult);
% Define fontsize for figures
fontsize_xylabel = 20 * t_mult;
fontsize_axis = 18 * t_mult;
fontsize_legend = 18 * t_mult;
% fontsize_text = 14 * t_mult;
% fontsize_table = 16 * t_mult;
% fontsize_titlepage = 30 * t_mult;

switch K
    case 1  % -> all data used for training
        plot(y,y_e,'b.','markersize',12); hold on; axis square
        plot(y_plot,y_plot,'k-','color',[0.5 0.5 0.5],'linewidth',2);
        axis tight;
        % Add labels
        xlabel('Actual y','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
        ylabel('Simulated y','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
        set(gca,'fontsize',fontsize_axis,'tickdir','out');
        % Add least squares line
        P = polyfit(y,y_e,1);
        plot(y,P(1)*y + P(2),'k-.','linewidth',2);
        % Add legend
        legend('{\color{red} Training}','{\color{gray} 1:1 Line}', ...
            '\color{black} Regression]', ...
            'location','SouthEast','box','off', ...
            'fontsize',fontsize_legend);
        % Add text
        RMSE_train = sqrt(1/N * sum((y - y_e).^2));
        % Now print results/statistics
        a = axis;
        text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.92*(a(4)-a(3)), ...
            strcat('Number of training data points:',{' '},...
            num2str(N)),'fontsize',fontsize_axis);
        text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.85*(a(4)-a(3)), ...
            strcat('RMSE for training data set:',{' '},...
            num2str(floor(T*RMSE_train)/T)),'fontsize',fontsize_axis);
    otherwise   % --> training and evaluation data set
        subplot(1,2,1),plot(y(id(1:R)),y_e(id(1:R)),'r+', ...
            'markersize',12,'linewidth',2); hold on; axis square
        plot(y_plot,y_plot,'k-','color',[0.5 0.5 0.5],'linewidth',2);
        axis tight;
        % Add labels
        xlabel('Actual $y$','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
        ylabel('Simulated $y$','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
        set(gca,'fontsize',fontsize_axis,'tickdir','out');
        % Add least squares line
        P = polyfit(y(id(1:R)),y_e(id(1:R)),1);
        plot(y(id(1:R)),P(1)*y(id(1:R))+P(2),'k-.','linewidth',2);
        % Add legend
        legend('{\color{red} Training}','{\color{gray} 1:1 Line}', ...
            '{\color{black} Regression}', ...
            'location','SouthEast','box','off', ...
            'fontsize',fontsize_legend);
        % Add text
        RMSE_train = sqrt(1/R * sum((y(id(1:R)) - y_e(id(1:R))).^2));
        % Now print results/statistics
        a = axis;
        text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.92*(a(4)-a(3)), ...
            strcat('Number of training data points:',{' '},...
            num2str(R)),'fontsize',fontsize_axis);
        text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.85*(a(4)-a(3)), ...
            strcat('RMSE for training data set:',{' '},...
            num2str(floor(T*RMSE_train)/T)),'fontsize',fontsize_axis);

        % Plot the evaluation data points
        id_all = ones(N,1); id_all(id) = 0; id_e = find(id_all>0);
        subplot(1,2,2),plot(y(id_e),y_e(id_e),'b.', ...
            'markersize',12); hold on; axis square
        plot(y_plot,y_plot,'k-','color',[0.5 0.5 0.5],'linewidth',2);
        axis tight;
        % Add labels
        xlabel('Actual $y$','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
        ylabel('Simulated $y$','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
        set(gca,'fontsize',fontsize_axis,'tickdir','out');
        % Add least squares line
        P = polyfit(y(id_e),y_e(id_e),1);
        plot(y(id_e),P(1)*y(id_e) + P(2),'k-.','linewidth',2);
        legend('{\color{blue} Evaluation}','{\color{gray} 1:1 Line}', ...
            '{\color{black} Regression}', ...
            'location','SouthEast','box','off', ...
            'fontsize',fontsize_legend);
        % Add text
        RMSE_eval = sqrt(1/(N-R) * sum((y(id_e) - y_e(id_e)).^2));
        % Now print results/statistics
        a = axis;
        text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.92*(a(4)-a(3)), ...
            strcat('Number of evaluation data points:',{' '},...
            num2str(N-R)),'fontsize',fontsize_axis);
        text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.85*(a(4)-a(3)), ...
            strcat('RMSE for evaluation data set:',{' '},...
            num2str(floor(T*RMSE_eval)/T)),'fontsize',fontsize_axis);
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
