%% Global Sensitivity Analysis of van Genuchton Water Retention Function %%
% In this script, we consider correlated parameters
clc;clear all;close all;

%% 1. Define Covariance Matrix
% Reference: Gomes et al (2017), The role of uncertainty in bedrock depth and hydraulic properties on the stability of a
% variably-saturated slope. Computers and Geotechnics, 88, 222–241. https://doi.org/10.1016/j.compgeo.2017.03.016
meanX = [0.4337 0.0577 0.0053 1.6391];
st = [0.0159 0.0153 0.001 0.055];

sigma = st.*ones(4).*st';              

C = eye(4); 
C(1,2) = 0.09; C(2,1) = 0.09;       % cov(theta_s, theta_r)
C(1,3) = 0.24; C(3,1) = 0.24;       % cov(theta_s, alpha)
C(1,4) = -0.04; C(4,1) = -0.04;     % cov(theta_s, n)
C(2,3) = 0.15; C(3,2) = 0.15;       % cov(theta_r, alpha)
C(2,4) = -0.58; C(4,2) = -0.58;     % cov(theta_r, n)
C(3,4) = -0.86; C(4,3) = -0.86;     % cov(alpha,n)
C = sigma.*C;              

%% 2. Global Sensitivity Analysis-HDMR_EXT
nh = 200;
h = -logspace(-1,6,nh);              % Model Input: Soil Water Pressure Head
N = 5000;                            % Number of Samples
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
V_sum = zeros(nh,1);
V_std = zeros(nh,d);

wt = waitbar(0,'HDMRext calculating');
for j = 1:nh
    Y = theta(:,j);
    % Random error
%     sigma2 = var(Y)/100;
    % Add random error
%     Y = Y + normrnd(0,sqrt(sigma2),N,1);
    options = struct('graphics',0,'basis',1,'maxorder',3,'maxiter',100,'bf1',1,'m',2,...
        'K',1,'R',N/2,'method',1,'alfa',0.99,'lambda',0.10,'vartol',1e-3,'refit',1);
    % Now run the HDMR_EXT toolbox
    [S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);
    % Extract the sensitivity indices
    for k = 1:n_ind
        S_a(j,k) = str2num(cell2mat(S(k+1,3)));
        S_b(j,k) = str2num(cell2mat(S(k+1,5)));
    end
    for k = 1:d
        S_t(j,k) = str2num(cell2mat(S(k+1,9)));
    end
    waitbar(j/nh,wt)
end

%% Make the plot
figure('name','figure 1','units','normalized','outerposition',[0.23 0.29 0.55 0.6])
parameters = ["$\theta_{\rm s}$","$\theta_{\rm r}$","$\alpha_{\rm vg}$","$n_{\rm vg}$"];
letter = 'a':'z';

for i = 1:d
    subplot(2,2,i)
    hold on
    p3 = semilogx(abs(h),S_t(:,i),'Color',[228, 98, 23]/255,'LineWidth',2.0);
    p1 = semilogx(abs(h),S_a(:,i),'--','Color',[106, 150, 57]/255,'LineWidth',2.5);
    p2 = semilogx(abs(h),S_b(:,i),'--','color',[57, 119, 219]/255,'LineWidth',2.5);
    p4 = semilogx(abs(h),S_t(:,i) - S_b(:,i) - S_a(:,i),'k-.','LineWidth',1.5);
    set(gca, 'XScale', 'log');
    % Set up the ticks
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'TickDir','Out')
    set(gca,'XMinorTick','on','TickLength',[0.02, 0.04])
    set(gca,'YMinorTick','on')
    set(gca,'LineWidth',1.0)
    set(gca,'fontsize',15)
    if i == 1 || i == 2
        set(gca,'Xticklabel',[]);
    end
    % Add x y labels
    if i == 3 || i == 4
        xlabel('$-h\, ({\rm cm})$','Interpreter','Latex','fontname','times','FontSize',18)
    end
    if i == 1 || i == 3
        ylabel('Sensitivity indices','Interpreter','Latex','fontname','times','FontSize',16)
    end
    ytickformat('%.1f')
    if i == 1
        bb = legend([p1,p2,p3,p4],' ',' ',' ',' ','Interpreter','Latex','fontname','times','fontsize',17);
        bb.Box = 'off';
        % add colored legend
        text(10^(0.10),1.05-0.19,'$S^{\rm a}_{i}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','Color',[106, 150, 57]/255)
        text(10^(0.10),0.88-0.19,'$S^{\rm b}_{i}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color',[57, 119, 219]/255)
        text(10^(0.10),0.71-0.19,'$S_{i}^{\rm T}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color',[228, 98, 23]/255)
        text(10^(0.10),0.55-0.19,'$\sum S_{ij} + \sum S_{ijk}$','Interpreter',...
        'Latex','fontname','times','fontsize',18,'horizontalAlignment',...
        'left','color','k')
    end
    % Set graph position
    if i <=2
        set(gca,'OuterPosition',[(i-1)*0.51 0.53 0.48 0.46])
    else
        set(gca,'OuterPosition',[(i-3)*0.51 0.01 0.48 0.51])
    end
    switch i
        case 1
            ylim([-0.2,1.2])
        case 2
            ylim([-0.2,1.2])
        case 3
            ylim([-1,1.5])
%             set(gca,'Ytick',-0.8:0.4:0.8)
        case 4
            ylim([-0.51,0.6])
    end
    ylmt = get(gca,'YLim');
%     xlmt = get(gca,'XLim');
    sublabel = ['(',letter(i),') ',char(parameters(i))];
    text(10^(-0.8),ylmt(1)+(ylmt(2) - ylmt(1))*0.92,sublabel,'Interpreter',...
        'Latex','fontsize',18,'horizontalAlignment','left')
    xlim([1e-1,1e6])

end
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% 
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters figure2
