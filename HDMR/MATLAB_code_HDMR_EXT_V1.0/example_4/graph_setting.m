function graph_setting(x_lab,y_lab,FigPar)
box off
set(gca,'TickLabelInterpreter','Latex')
set(gca,'TickDir','Out')
set(gca,'fontsize',FigPar.gca_ftsize,'fontname','times')
set(gca,'YMinorTick','on','XMinorTick','on')
% set(gca,'XMinorTick','on')
set(gca,'LineWidth',0.8)
if isempty(x_lab) == 0
    xlabel(x_lab,'interpreter','latex','fontsize',FigPar.axislabel_size);
end

if isempty(y_lab) == 0
    ylabel(y_lab,'interpreter','latex','fontsize',FigPar.axislabel_size);
end

ax = gca;
ax.TickLength = FigPar.ax_ticklength;