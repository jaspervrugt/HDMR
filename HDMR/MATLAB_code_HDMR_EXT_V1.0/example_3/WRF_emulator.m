% In this script, we show the emulation capacity of HDMR_EXT
clear;clc

%% RGB color
color_RGB = [245, 179, 66; 164, 245, 66; 69, 245, 66; 66, 245, 197;...
    66, 197, 245; 66, 117, 245; 96, 66, 245; 173, 66, 245;...
    245, 66, 209; 245, 66, 108; 156, 129, 135; 79, 58, 62];

%%
% Load the VG parameters
VG_pars=xlsread('MVG_pars.xlsx','Carsel, Parish','B3:F14');
% Define soil texture
soil_texture=["Sand","Loamy sand","Sandy loam","Sandy clay loam","Loam","Sandy clay",...
    "Silt loam","Silt","Clay loam","Silty clay loam","Silt clay","Clay"];
% Parameter boundary
theta_s = [0.2,0.5];                % Saturated Volumetric Moisture Content [m3/m3]
theta_r = [0,0.1];                  % Residual Volumetric Moisture Content  [m3/m3]
alpha = [0.001,0.3];                 % Reciprocal of Air Entry Value [1/m]
n = [1.1, 4.0];                     % A Parameter which is Proportional to the Slope of Water Retention Function
par_bound = [theta_s; theta_r; alpha; n]';

% Number of samples
N = 5000;
d = size(VG_pars,2);
nh = 100;

st = [0.0159 0.0153 0.001 0.055];
sigma = st.*ones(4).*st';
C = eye(4); 
C = sigma.*C;

%% HDMR_EXT calibration

% Loop starting from the first soil
for i = 1:size(VG_pars,1)
    % Derive the input
    x = [VG_pars(i,2),VG_pars(i,1),VG_pars(i,3),1/(1-VG_pars(i,4))];
    % Random sampling
    X = mvnrnd(x,C,N);
    % Normalize the values
    X = (X - min(X))./(max(X) - min(X)) .* (par_bound(2,:)-par_bound(1,:)) + par_bound(1,:);
    N = size(X,1);
    X(1,:) = x;
    
    h = -logspace(-1,4,nh);
    % Preallocation
    theta = zeros(N,nh);
    for j = 1:N
        theta(j,:) = VG(X(j,:),h);
    end
    
    for j = 1:size(h,2)
        Y_m(i,j) = theta(1,j);
        Y = theta(:,j);
        % Random error
%         sigma2 = var(Y)/100;
        % Add random error
%         Y = Y + normrnd(0,sqrt(sigma2),N,1);
        options = struct('graphics',0,'maxorder',3,'maxiter',100,'bf1',1,'m',3,...
            'K',1,'R',N,'method',1,'alfa',0.99,'lambda',0.10,'vartol',1e-3,'refit',1);
        % Now run the HDMR_EXT toolbox
        [S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);
        vg_e(i,j) = Em.Y_e(1,1);
    end
end

%%
% Plot model versus emulator
fig = figure('name','Sensitivity Indices wrt h','units','normalized','outerposition',[0.2 0.2 0.39 0.56]);
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
lmt = 0:0.1:0.6; a = -1:44;
letter = 'a':'z';
y_tick = 10.^(a);

for i = 1:size(VG_pars,1)
    % Determine the upper and lower limit of theta
    tr = VG_pars(i,1);  ts = VG_pars(i,2);
    lw_bd = lmt(lmt < tr); lw_bd = lw_bd(1);
    up_bd = lmt(lmt > ts); up_bd = up_bd(1);
%     subplot(1,3,i)
    hold on
    pm(i) = semilogy(Y_m(i,:),abs(h),'o','markerfacecolor',[90 90 90]/255,'MarkerEdgeColor',[90 90 90]/255,'linewidth',0.5,'markersize',2);
    pd = semilogy(vg_e(i,:),abs(h),'-','color',[90 90 90]/255,'LineWidth',1.0);
   
    set(gca, 'YScale', 'log');
    % Determine the upper limit of h
    h_end = abs(h(end));
    if h_end < 1e+3
        set(gca,'YTick',y_tick(1:5))
        ylim([y_tick(1), y_tick(5)])
    elseif h_end > 1e+3 && h_end < 1e+4
        set(gca,'YTick',y_tick(1:2:7))
        ylim([y_tick(1), y_tick(6)])
    elseif h_end > 1e+4 && h_end < 1e+5
        set(gca,'YTick',y_tick(1:2:7))
        ylim([y_tick(1), y_tick(7)])
    elseif h_end > 1e+5 && h_end < 1e+8
        set(gca,'YTick',y_tick(1:3:10))
        ylim([y_tick(1), y_tick(10)])
    elseif h_end > 1e+8 && h_end < 1e+10
        set(gca,'YTick',y_tick(1:4:13))
        ylim([y_tick(1), y_tick(13)])
    elseif h_end > 1e+10 && h_end < 1e+11
        set(gca,'YTick',y_tick(1:7:15))
        ylim([y_tick(1), y_tick(13)])
    elseif h_end > 1e+11 && h_end < 1e+12
        set(gca,'YTick',y_tick(1:7:15))
        ylim([y_tick(1), y_tick(15)])
    elseif h_end > 1e+13 && h_end < 1e+30
        set(gca,'YTick',y_tick(1:4:17))
        ylim([y_tick(1), y_tick(17)])
    elseif h_end > 1e+30
        set(gca,'YTick',y_tick(1:10:41))
        ylim([y_tick(1), y_tick(41)])
        q = -1:2:39; 
        ax = gca; ax.YAxis.MinorTickValues = 10.^q;
    end
    
    set(gca,'YtickLabel',{'$10^{-1}$','$10^{0}\,$ ','$10^{1}\,$ ', '$10^{2}\,$ ', '$10^{3}\,$ ', '$10^{4}\,$ '})

    % Set up the ticks
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'TickDir','Out')
    set(gca,'XMinorTick','on','TickLength',[0.01, 0.02])
    set(gca,'YMinorTick','on')
    set(gca,'LineWidth',1.0)
    set(gca,'fontsize',14)
    xtickformat('%.2f')
     
    % Add x,y label
    xlabel('Volumetric soil moisture content, $\theta \,({\rm cm^{3}/cm^{3}})$','Interpreter','Latex','fontsize',16)
    ylabel('$-$ Soil water pressure head, $-h \,({\rm cm})$','Interpreter','Latex','fontsize',16)
end


%% Add the text
letters = 'a':'l';

text(0.071447366632913,22.368,'a','Interpreter','latex','fontsize',14)
text(0.10,50,'b','Interpreter','latex','fontsize',14)
text(0.125394738034198,111.099,'c','Interpreter','latex','fontsize',14)
text(0.268552633963133,31.25,'d','Interpreter','latex','fontsize',14)
text(0.285526327710403,51.03,'e','Interpreter','latex','fontsize',14)
text(0.309473677057969,71.23,'f','Interpreter','latex','fontsize',14)
text(0.242105263157895,261.48,'g','Interpreter','latex','fontsize',14)
text(0.209736831094089,716.65,'h','Interpreter','latex','fontsize',14)
text(0.170973689932572,3151.29,'i','Interpreter','latex','fontsize',14)
text(0.229105266614964,3151.29,'j','Interpreter','latex','fontsize',14)
text(0.280131570602718,3151.29,'k','Interpreter','latex','fontsize',14)
text(0.3101,3151.29,'l','Interpreter','latex','fontsize',14)
k = 1;
for i = 1:3
    for j = 1:4
           lgd = ['(',letters(k),') '];
           text(0.01 + (i-1) * 0.12, 10^(0.5 - (j-1) * 0.4), lgd,'color','k', 'Interpreter','latex','fontsize',12)
           text(0.03 + (i-1) * 0.12, 10^(0.5 - (j-1) * 0.4), soil_texture(k),'color','k', 'Interpreter','latex','fontsize',12)
           k = k + 1;
    end
end
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% 
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters HDMRext_emulator_WRF
