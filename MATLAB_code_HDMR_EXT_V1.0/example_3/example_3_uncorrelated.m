%% Global Sensitivity Analysis on Van Genuchton Water Retention Function %%
% In this script, we consider uncorrelated parameters
clc;clear all;
% Global Sensitivity Analysis on Van Genuchton Water Retention Function

%% 1. Define Covariance Matrix
% Define the mean value: Gomes et al (2017)
meanX = [0.4337 0.0577 0.0053 1.6391];
st = [0.0159 0.0153 0.001 0.055];

sigma = st.*ones(4).*st';

% Define covariance matrix: 
C = eye(4); 
C = sigma.*C;

%% 2. Global Sensitivity Analysis-HDMR_EXT
nh = 100;
h = -logspace(-1,6,nh);              % Model Input: Soil Water Pressure Head
N = 5000;                            % Number of Samples

% Generate the samples
X = mvnrnd(meanX,C,N);

% How many parameters?
d = size(X,2); 
for i = 1:N
    theta(i,:) = VG(X(i,:),h);
end
% Preallocation
n_ind = d + nchoosek(d,2) + nchoosek(d,3);
S_a = zeros(nh,n_ind);
S_b = zeros(nh,n_ind);
S_t = zeros(nh,d);
V_std = zeros(nh,d);

for j = 1:size(h,2)
    Y = theta(:,j);
    % Random error
%     sigma2 = var(Y)/100;
    % Add random error
%     Y = Y + normrnd(0,sqrt(sigma2),N,1);
    options = struct('graphics',0,'basis',1,'maxorder',3,'maxiter',100,'bf1',1,'m',3,...
        'K',1,'R',N/2,'method',1,'alfa',0.99,'lambda',0.10,'vartol',1e-3,'refit',1);
    % Now run the HDMR_EXT toolbox
    [S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);
    vg_e(i,j) = Em.Y_e(1,1);
    % Extract the sensitivity indices
    for k = 1:n_ind
        S_a(j,k) = str2num(cell2mat(S(k+1,3)));
        S_b(j,k) = str2num(cell2mat(S(k+1,5)));
    end
    for k = 1:d
        S_t(j,k) = str2num(cell2mat(S(k+1,9)));
%         V_std(j,k) = str2num(cell2mat(S(k+1,12)));
    end
end

%% Jacobian determination
x = meanX;         
theta_old = VG(x,h)';
x_old = x;
% Derive the numeric partial derivative
for j = 1:4
   x = x_old;
   delta_x = 0.001*x(j); 
   x(j) = x(j) + delta_x;
   theta_new = VG(x,h)'; 
   J(:,j) = (theta_new - theta_old)/delta_x;
end

%% Make the plot
figure('name','Jacobian wrt h','units','normalized','outerposition',[0.23 0.29 0.57 0.6])
parameters = ["$\theta_{\rm s}$","$\theta_{\rm r}$","$\alpha_{\rm vg}$","$n_{\rm vg}$"];
letter = 'a':'z';

for i = 1:d
    subplot(2,2,i)
    hold on
    % The left axis -- global sensitivity estimates
    yyaxis left
    if i == 1 || i == 3
        ylabel('Sensitivity indices','Interpreter','Latex','fontname','times','FontSize',14)
    end
    
    p3 = semilogx(abs(h),S_t(:,i),'Color',[228, 98, 23]/255,'LineWidth',2.0);
    p1 = semilogx(abs(h),S_a(:,i),'--','Color',[106, 150, 57]/255,'LineWidth',2.5);
    p2 = semilogx(abs(h),S_b(:,i),'--','color',[57, 119, 219]/255,'LineWidth',2.5);
    set(gca, 'XScale', 'log');
    ax = gca; % current axes
    ax.FontSize = 12;
    set(gca,'TickDir','Out')
    set(gca,'XMinorTick','on','TickLength',[0.02, 0.04])
    set(gca,'YMinorTick','on')
    set(gca,'LineWidth',1.0)
    set(gca,'ycolor','k')
    ytickformat('%.1f')
    switch i
        case 1
            ylim([-0.2,1.2])
        case 2
            ylim([-0.2,1.2])
        case 3
            ylim([-0.02,0.8])
        case 4
            ylim([-0.02,0.4])
    end
    % The right axis -- local sensitivity estimates
    yyaxis right
    p4 = semilogx(abs(h),abs(J(:,i)),'-','color','k','LineWidth',1.5);
    if i == 2 || i == 4
        ylabel('Analytical local sensitivity','Interpreter','Latex','fontname','times','FontSize',14)
    end
    ytickformat('%.1f')
    set(gca, 'XScale', 'log');
    % Set up the ticks
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'TickDir','Out')
    set(gca,'XMinorTick','on','TickLength',[0.02, 0.04])
    set(gca,'YMinorTick','on')
    set(gca,'LineWidth',1.0)
    set(gca,'fontsize',15)
    set(gca,'ycolor','k')
    if i == 1
        bb = legend([p1,p2,p3,p4],' ',' ',' ',' ','Interpreter','Latex','fontname','times','fontsize',18);
        bb.Box = 'off';
        % add colored legend
        text(10^(0.2),1.22-0.5,'$S^{\rm a}_{i}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color',[106, 150, 57]/255)
        text(10^(0.2),1.04-0.5,'$S^{\rm b}_{i}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color',[57, 119, 219]/255)
        text(10^(0.2),0.86-0.5,'$S_{i}^{\rm T}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color',[228, 98, 23]/255)
        text(10^(0.2),0.69-0.5,'$|\partial \theta_{j}/ \partial x_{i}|$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color','k')
    end
    % Set graph position
    if i <=2
        set(gca,'OuterPosition',[(i-1)*0.51 0.53 0.48 0.46])
    else
        set(gca,'OuterPosition',[(i-3)*0.51 0.01 0.48 0.51])
    end
    if i == 3 || i == 4
        xlabel('$-h\, ({\rm cm})$','Interpreter','Latex','fontname','times','FontSize',18)
    end
    if i == 1 || i == 2
        set(gca,'Xticklabel',[]);
    end
    
    switch i
        case 1
            ylim([-0.2,1.2])
        case 2
            ylim([-0.2,1.2])
        case 3
            ylim([-0.1,20])
%             set(gca,'Ytick',-0.8:0.4:0.8)
            ytickformat('%.0f')
        case 4
            ylim([-0.001,0.3])
    end
    ylmt = get(gca,'YLim');
%     xlmt = get(gca,'XLim');
    sublabel = ['(',letter(i),') ',char(parameters(i))];
    text(10^(-0.8),ylmt(1)+(ylmt(2) - ylmt(1))*0.92,sublabel,'Interpreter','Latex','fontname','times','fontsize',18,'horizontalAlignment','left')
    xlim([1e-1,1e6])
end

% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% 
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters WRC_HDMR-EXT(Uncorrelated)_revised


