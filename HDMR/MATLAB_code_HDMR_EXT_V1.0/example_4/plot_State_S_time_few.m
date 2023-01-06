function plot_State_S_time_few(S_indices,S_sim,Y_MAP)

FigPar.gca_ftsize = 10;
FigPar.axislabel_size = 12;
FigPar.ax_ticklength = [0.005 0.01];
%% Load measured data
load('daily_data_maurer_1980_2008.mat')
t_daily = datetime(daily_data(:,3),daily_data(:,2),daily_data(:,1)); 
P_daily = daily_data(:,4);    % Precipitation in mm/d
Ep_daily = daily_data(:,5);   % Potential evaporation in mm/d
Q_daily = daily_data(:,6);    % Discharge in mm/d

% Calibration period
t1 = datetime(1999,10,01);
t2 = datetime(2008,09,30);
t = t1:t2;
idx_t1 = find(t_daily==t1);
idx_t2 = find(t_daily==t2);
letter = 'a':'z';
idx1 = 367; idx2 = 732;

%% Derive the state variables
Int = S_sim(2:end,1);          % Incerception storage
S_uz = S_sim(2:end,5);         % Unsaturated zone storage
S_fr = S_sim(2:end,6);         % Fast reservoir storage
S_sr = S_sim(2:end,12);        % Slow reservoir storage

% %% Or two day average?
% Int = (S_sim(1:end-1,1) + S_sim(2:end,1))/2;          % Incerception storage
% S_uz = (S_sim(1:end-1,5) + S_sim(2:end,5))/2;         % Unsaturated zone storage
% S_fr = (S_sim(1:end-1,6) + S_sim(2:end,6))/2;         % Fast reservoir storage
% S_sr = (S_sim(1:end-1,12) + S_sim(2:end,12))/2;        % Slow reservoir storage
% 
% %% Or two day difference?
% Int = (S_sim(1:end-1,1) - S_sim(2:end,1))/2;          % Incerception storage
% S_uz = (S_sim(1:end-1,5) - S_sim(2:end,5))/2;         % Unsaturated zone storage
% S_fr = (S_sim(1:end-1,6) - S_sim(2:end,6))/2;         % Fast reservoir storage
% S_sr = (S_sim(1:end-1,12) - S_sim(2:end,12))/2;        % Slow reservoir storage

%% 
parameters = ["$R_{\rm u,max}$","$Q_{\rm s,max}$","$\alpha_{\rm e}$","$\alpha_{\rm f}$","$K_{\rm f}$"];

figure('name','Cluster','units','normalized','outerposition',[0.3 0.02 0.40 0.75])

par_idx = [2 3 4 5 6];

num_row = 8;

% Precipitation
for i = 1:num_row
    
    switch i
        case 1 % precipitation
            ax1 = axes('Position',[0.08 0.87 - (i - 1)*0.12, 0.82,0.1]);
            plot(t,Y_MAP,'color',[228, 98, 23]/255,'LineWidth',1.5)
            xlim([t(idx1) t(idx2)])
            graph_setting([],'$Q_{\rm t}$ (mm/d)',FigPar)
            set(gca,'xticklabel',[])
            ylim([-5,45])
            hAx(1)=gca;
            hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
            hold(hAx(2),'on')
            bar(t_daily,-P_daily,'facecolor','b')
            graph_setting([],'$P_{t}$ (mm/d)',FigPar)
            set(gca,'yticklabel',{'100','50','0'},'xtick',[])
            ylim([-100,0])
        case 2
            ax1 = axes('Position',[0.08 0.87 - (i - 1)*0.12, 0.82,0.1]);
            plot(t,S_uz,'color',[228, 98, 23]/255,'LineWidth',1.5)
            graph_setting([],'$R_{\rm u}$ (mm)',FigPar)
             ylim([-40,440])
             set(gca,'Ytick',0:200:400,'YtickLabel',{'$0$','$200$','$\,\,\,\,400$'})
        case 3
            ax1 = axes('Position',[0.08 0.87 - (i - 1)*0.12, 0.82,0.1]);
            plot(t,S_fr,'color',[228, 98, 23]/255,'LineWidth',1.5)
            graph_setting([],'$R_{\rm f}$ (mm)',FigPar)
            ylim([-10,110])
        case {4,5,6,7,8}
            hold on
            % Sb
            ax2 = axes('Position',[0.08 0.87 - (i - 1)*0.12, 0.82,0.1]);
            plot(t,S_indices.Sb(:,par_idx(i-3)),':','Color',[57, 119, 219]/255,'LineWidth',1.5)
            set(ax2,'ycolor',[57, 119, 219]/255)
            ax2.YAxisLocation = 'right';
            ytickformat('%.1f')
            switch i
                case 4
                    ylim([-0.8,0.4])
                    set(ax2,'Ytick',-0.6:0.4:0.2,'YtickLabel',{'$-0.6$','$-0.2$','$\,\,\,\,0.2$'})
                case 5
                    ylim([-0.4,0.4])
                    set(ax2,'Ytick',-0.3:0.3:0.3,'YtickLabel',{'$-0.3$','$\,\,\,\,0.0$','$\,\,\,\,0.3$'})
                case 6
                    ylim([-2.0,1.0])
                    set(ax2,'Ytick',-1.4:1.0:1.4,'YtickLabel',{'$-1.4$','$-0.4$','$\,\,\,\,0.4$','$\,\,\,\,1.4$'})
                case 7
                    ylim([-2.0,1.0])
                    set(ax2,'Ytick',-1.4:1.0:1.4,'YtickLabel',{'$-1.4$','$-0.4$','$\,\,\,\,0.4$','$\,\,\,\,1.4$'})
                case 8
                    ylim([-0.3,0.08])
                    set(ax2,'Ytick',-0.24:0.12:0.00,'YtickLabel',{'$-0.24$','$-0.12$','$\,\,\,\,0.00$'})
            end
            graph_label = ['$S^{\rm b}_{',num2str(par_idx(i-3)),'}$'];
            graph_setting([],graph_label,FigPar)
            xlim([t(idx1) t(idx2)])
            % Sa
            ax1 = axes('Position',[0.08 0.87 - (i - 1)*0.12, 0.82,0.1]);
            plot(t,S_indices.Sa(:,par_idx(i-3)),'Color',[106, 150, 57]/255,'LineWidth',1.5)
            set(ax1,'ycolor',[106, 150, 57]/255)
            graph_label = ['$S^{\rm a}_{',num2str(par_idx(i-3)),'}$'];
            graph_setting([],graph_label,FigPar)
            ytickformat('%.1f')
            switch i
                case {4,8}
                    ylim([-0.2,1])
                    set(ax1,'Ytick',0:0.4:0.8,'YtickLabel',{'$0.0$','$0.4$','$0.8$'})
                case 5
                    ylim([-0.3,1.3])
                    set(ax1,'Ytick',-0.5:0.5:1.5,'YtickLabel',{'$-0.5$','$0.0$','$0.5$','$1.0$','$1.5$'})
                case 7
                    ylim([-0.5,2.6])
                case 6
                    ylim([-0.4,2.4])
            end
            xlim([t(idx1) t(idx2)])
            set(gca, 'Color', 'None')
            if i ~= num_row
                set(ax1,'XtickLabel',[])
                set(ax2,'XtickLabel',[])
            end
    end
    if i ~= num_row
        set(gca,'XtickLabel',[])
    end
    xlim([t(idx1) t(idx2)])
    % Add the figure label
    if i > 3
        fig_label = ['(',letter(i),') ',char(parameters(i-3))];
    else
        fig_label = ['(',letter(i),') '];
    end
    xlim([t(idx1) t(idx2)])
    xlmt = get(gca,'XLim');
    ylmt = get(gca,'YLim');
    % add a border
    if i == 2||i == 3
        hold on
        plot([xlmt(2) xlmt(2)], [ylmt(1) ylmt(2)],'k','linewidth',0.8)
        plot([xlmt(1) xlmt(2)], [ylmt(2) ylmt(2)],'k','linewidth',0.8)
    elseif i ~= 1
        hold on
        plot([xlmt(1) xlmt(2)], [ylmt(2) ylmt(2)],'k','linewidth',0.8)
    end
        
    text(xlmt(1)+(xlmt(2) - xlmt(1))*0.02,ylmt(1)+(ylmt(2) - ylmt(1))*0.88,fig_label,'Interpreter','Latex',...
        'fontsize',FigPar.axislabel_size-1,'horizontalAlignment','left')
end
    