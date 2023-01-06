function plot_S_prior_posterior(S_indices_prior, S_indices_unc, S_indices)

%% Load measured data
t1 = datetime(1999,10,01);
t2 = datetime(2008,09,30);
t = t1:t2;
idx1 = 367; %idx1 = 1462;
idx2 = 732; %idx2 = 1827;

include_S = 1;

%% Set up graph information
FigPar.axislabel_size = 14;
legend_size = 14;
FigPar.gca_ftsize = 12;
Line_Width = 1.5;
FigPar.ax_ticklength = [0.005 0.010];
date_tick_idx = [398,459,518,579,640,702];
X_tick = t(date_tick_idx);

%% Plot model simulation with 95% confidence interval
% Model parameters
parameters = ["$S_{\rm max}$","$\alpha_{\rm E}$",...
    "$\alpha_{\rm F}$","$K_{\rm F}$","$K_{\rm S}$"];
parameters = ["$I_{\rm max}$","$R_{\rm u,max}$","$Q_{\rm s, max}$","$\alpha_{\rm e}$",...
    "$\alpha_{\rm f}$","$K_{\rm f}$","$K_{\rm s}$"];

letter = 'a':'z';

figure('name','Parameter sensitivity as function time','units','normalized','outerposition',[0.2 0.10 0.4 0.85])


%% Plot sensitivity indices as a function of time

for i = 1:7
    ax1 = axes('Position',[0.09,0.86 - (i-1)*0.138,0.86,0.11]);
    hold on
    plot(t,S_indices_prior.St(:,i),':','Color',[144, 66, 245]/255,'LineWidth',Line_Width)
    plot(t,S_indices_unc.St(:,i),'Color',[247, 153, 12]/255,'LineWidth',Line_Width)
    if include_S == 1
%         plot(t,S_indices.Sa(:,i)+S_indices.Sb(:,i),'-.','Color',[106, 150, 57]/255,'LineWidth',Line_Width)
        plot(t,S_indices.St(:,i),'-.','Color',[106, 150, 57]/255,'LineWidth',Line_Width)
    end
    graph_label = ['$S_{',num2str(i),'}$'];
    if include_S == 1
        graph_setting([],graph_label,FigPar)
    else
        graph_setting([],graph_label,FigPar)
    end
    ytickformat('%.1f')
    switch i
        case 1
            ylim([-0.1,1.0])
        case 2
            ylim([-0.3,1.0])
        case {6,7}
            ylim([-0.1,1])
        case 3
            ylim([-0.1,1.0])
        case 4
            ylim([-0.4,1.5])
        case 5
            ylim([-0.3,1])
    end
    if i == 1
        if include_S == 1
%             plot(t(630:645),1.46*(ones(1,16)),':','color',[144, 66, 245]/255,'LineWidth',Line_Width);
%             plot(t(630:645),1.10*(ones(1,16)),'-','color',[247, 153, 12]/255,'linewidth',1.5);
%             plot(t(630:645),0.74*(ones(1,16)),'-.','color',[106, 150, 57]/255,'linewidth',1.5);
            legend('','','','fontsize',12,'interpreter','latex')
            legend box off
            text(t(650),1.50-0.3,'Prior hypercube','interpreter','latex','fontsize',12,'color',[144, 66, 245]/255)
            text(t(650),1.28-0.3,'Posterior hypercube','interpreter','latex','fontsize',12,'color',[247, 153, 12]/255)
            text(t(650),1.06-0.3,'Posterior distribution','interpreter','latex','fontsize',12,'color',[106, 150, 57]/255)
        end
    end
    
    if i ~= 7
        set(gca,'XtickLabel',[])
    end

    % Add the figure label
    fig_label = ['(',letter(i),') ',char(parameters(i))];
    xlim([t(idx1) t(idx2)])
    xlmt = get(gca,'XLim');
    ylmt = get(gca,'YLim');
    text(xlmt(1)+(xlmt(2) - xlmt(1))*0.01,ylmt(1)+(ylmt(2) - ylmt(1))*0.88,fig_label,'Interpreter','Latex',...
        'fontsize',FigPar.axislabel_size-1,'horizontalAlignment','left')
end