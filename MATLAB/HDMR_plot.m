function HDMR_plot(SA_sig,Fx,y,y_e,select,p0,iter,id,R,K,refit,maxorder)
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
% Plot the results of HDMR toolbox                                        %
%                                                                         %
%  SYNOPSIS                                                               %
%   HDMR_plot(SA_sig,Fx,y,y_e,select,p0,iter,id,R,K,refit,maxorder)       %
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Print wait statement to the screen
fprintf('HDMR PLOTTING: PLEASE WAIT ...');
% Define name of program
n_program = 'HDMR';
% Define name of figures file
file_name = [n_program,'_figures.pdf'];
% Define names of tables
table1_name = ['Table_1_',n_program,'.pdf']; 
table2_name = ['Table_2_',n_program,'.pdf'];
% Determine the entries without pdf
id_pdf = 1:(strfind(table1_name,'pdf') - 2);

% import mlreportgen.dom.*;
% doc = Document('HDMR_figures','pdf'); open(doc)
% Determine screen size
scr_z = get(0,'ScreenSize'); 
% Multiplier, x and y: axis
x_mult = scr_z(3)/1920; y_mult = scr_z(4)/1080;
% Multiplier, text
t_mult = min(x_mult,y_mult);
% Define fontsize for figures
fontsize_xylabel = 18 * t_mult;
fontsize_axis = 16 * t_mult;
fontsize_legend = 14 * t_mult;
fontsize_text = 14 * t_mult;
fontsize_table = 16 * t_mult;
fontsize_titlepage = 30 * t_mult;

% ----------------------------------------------------------------------- %
%                 Now plot empty figure for PDF file                      %
% ----------------------------------------------------------------------- %

figure('units','normalized','outerposition',[0 0 1 1],'name','first page');
plot([],[],'ro'); axis([0 1 0 1]); set(gcf,'color','w'); 
set(gca,'XColor','w','YColor','w');
%title('Visual results of HDMR toolbox','fontsize',20,'interpreter','latex');
text(0.3*x_mult,0.6*y_mult,'Visual results of HDMR toolbox','fontsize', ...
    fontsize_titlepage,'interpreter','latex');
% Now about Tables
text(0.3*x_mult,0.5*y_mult,'$\;\;$GUI Tables may not print well ', ...
    'fontsize', fontsize_titlepage,'interpreter','latex');

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                   Plot results of emulator of K trials                  %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% First plot scatter plot with results of each emulator
y_plot = linspace(min(y),max(y));
% plot emulators in subplots with r rows and c columns
row = 2; col = 5;
% Now loop over each emulator and make scatter plot
for k = 1 : K
    % If true: open new figure and reset nr of columns to one
    if (mod(k,row*col) == 1) || (k == 1)
        figure('units','normalized','outerposition',[0 0 1 1],...
            'name',['HDMR: Bivariate scatter plots actual and ' ...
            'emulated y values of training data set']); c = 1;
    end
    % determine appropriate index of figure
    if c <= col
        r = 1; plot_c = c;
    else
        r = 2; plot_c = c - col;
    end
    % select right Y indices
    id_R = id(1:R,k);
    % Define the location of the figure
    ax1 = axes('units','inches');
    % define new axis position
    axpos1 = [0.9 + (plot_c - 1) * 3.75 6 - (r-1) * 4.75 3 3]; 
    % scale in x and y
    axpos1([1 3]) = axpos1([1 3]) * x_mult;
    axpos1([2 4]) = axpos1([2 4]) * y_mult;
    % End scale in x and y
    set(ax1,'position',axpos1);
    % Now plot the training data
    plot(y(id_R),y_e(1:R,k),'s','MarkerFaceColor','r', ...
        'MarkerEdgeColor','r'); hold on; axis square; axis tight;
    % Now plot the 1:1 Line in gray
    plot(y_plot,y_plot,'k-','color',[0.5 0.5 0.5], ...
        'linewidth',2); axis tight;
    % Add labels
    xlabel('$y$','fontsize',fontsize_xylabel,'interpreter','latex');
    if mod(k,col) == 1
        ylabel('$y = f(\textbf{x})$','fontsize',fontsize_xylabel, ...
            'interpreter','latex');
    end
    set(gca,'fontsize',fontsize_axis);
    % add legend only in top-left plot
    if ( c == 1 )
        leg = legend({'Training','1:1 Line'},'location', ...
            'southeast','interpreter','latex');
        set(leg,'linewidth',3); set(leg,'fontsize',fontsize_legend); 
        legend boxoff;
    end
    % Calculate training statistics
    RMSE_train = sqrt(1/R * sum((y(id_R) - y_e(1:R,k)).^2));
    RMSE_1 = num2str(RMSE_train); RMSE_1 = RMSE_1(1:min(numel(RMSE_1),5));
    title(strcat('Bootstrap trial:',{' '},num2str(k)), ...
        'fontweight','bold','interpreter','latex');
    % Now print results/statistics
    a = axis;
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.92*(a(4)-a(3)), ...
        strcat('$\# terms.:',{' '},...
        num2str(sum(select(:,k))),'$'),'interpreter','latex', ...
        'fontsize',fontsize_text);
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.85*(a(4)-a(3)), ...
        strcat('$\# coef.:',{' '},...
        num2str(p0(k)),'$'),'interpreter','latex','fontsize',fontsize_text);
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.78*(a(4)-a(3)), ...
        strcat('$\# iter.:',{' '},...
        num2str(sum(iter(1:4,k))),'$'),'interpreter','latex', ...
        'fontsize',fontsize_text);
    text(a(1) + 0.05*(a(2)-a(1)),a(3) + 0.70*(a(4)-a(3)), ...
        strcat('${\rm RMSE}_{\rm T}:',{' '},RMSE_1,'$'),...
        'interpreter','latex','fontsize',fontsize_text);
    % update column number
    c = c + 1;
end
% Get figure handles
figHandles = flipud(findall(0,'Type','figure'));
for zz = 1:numel(figHandles)
    figure(figHandles(zz)); set(gcf,'color','w');
    switch zz
        case 1
            exportgraphics(figHandles(1),file_name);
        otherwise
            % Append to existing figure
            exportgraphics(figHandles(zz),file_name,'Append',true);
    end
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                   Plot number of required iterations                    %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

if K > 1
    % open figure
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.5],...
        'name',['Total number of iterations required for each ' ...
        'order as function of bootstrap trial']);
    % Define the location of the figure
    ax1 = axes('units','inches');
    % define new axis position
    axpos1 = [ 1.0 0.9 16.2 4 ]; 
    % scale in x and y
    axpos1([1 3]) = axpos1([1 3]) * x_mult;
    axpos1([2 4]) = axpos1([2 4]) * y_mult;
    % end scale in x and y
    set(ax1,'position',axpos1, ...
        'units','normalized');
    % now plot iteration versus # emulator
    plot(1:K,iter(1,1:K),'r','linewidth',2); hold on;
    h = plot(1:K,iter(1,1:K),'s','MarkerFaceColor','r', ...
        'MarkerEdgeColor','r');
    set(get(get(h,'Annotation'),'LegendInformation'), ...
        'IconDisplayStyle','off');
    print_lgd = 1;
    if ( maxorder > 1 )     % add second order
        plot(1:K,iter(2,1:K),'b','linewidth',2); hold on;
        h = plot(1:K,iter(2,1:K),'o','MarkerFaceColor','b', ...
            'MarkerEdgeColor','b');
        set(get(get(h,'Annotation'),'LegendInformation'), ...
            'IconDisplayStyle','off');
        print_lgd = [ 1 2 ];
    end
    if ( maxorder == 3 )    % add third order
        plot(1:K,iter(3,1:K),'g','linewidth',2); hold on;
        h = plot(1:K,iter(3,1:K),'d','MarkerFaceColor','g', ...
            'MarkerEdgeColor','g');
        set(get(get(h,'Annotation'),'LegendInformation'), ...
            'IconDisplayStyle','off');
        print_lgd = [ 1 2 3 ];
    end
    if ( refit == 1 )       % add refit
        plot(1:K,iter(4,1:K),'k','linewidth',2); hold on;
        h = plot(1:K,iter(4,1:K),'^','MarkerFaceColor','k', ...
            'MarkerEdgeColor','k');
        set(get(get(h,'Annotation'),'LegendInformation'), ...
            'IconDisplayStyle','off');
        print_lgd = [ print_lgd 4];
    end
    % determine axis
    a = axis; axis([0.5 K + 0.5 a(3) - 0.05*(a(4)-a(3)) ...
        a(4) + 0.05*(a(4)-a(3))]);
    % add title
    ylabel('Number of iterations','fontsize',fontsize_xylabel, ...
        'interpreter','latex');
    % Add labels
    xlabel('Emulator (bootstrap trial)','fontsize',fontsize_xylabel, ...
        'interpreter','latex');
    % add legend
    lgd = {'first order','second order','third order','refitting'};
    legend(lgd(print_lgd),'interpreter','latex');
    % fontsize
    set(gca,'fontsize',fontsize_axis); set(gcf,'color','w');
    % Add to PDF file
    exportgraphics(gca, file_name, 'Append', true);
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                    NEW: Plot Fx to table in figure                      %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Write table
% % import mlreportgen.dom.*;
% % t = Table(Fx(2:K+1,1:4));
% % t.Style = {RowHeight("1cm")};
% % t.Border = "solid";
% % t.BorderWidth = "1px";
% % t.ColSep = "solid";
% % t.ColSepWidth = "1";
% % t.RowSep = "solid";
% % t.RowSepWidth = "1";
% % t.TableEntriesStyle = [t.TableEntriesStyle {FontFamily("Arial"), ...
% %     Width("2cm"),Color("blue"),Bold}];
% % t.TableEntriesHAlign = "center"; t.TableEntriesVAlign = "middle";
% % doc = Document(file_name,'pdf');
% % append(doc,t); close(doc);

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                       Plot Fx to table in figure                        %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

col_names = Fx(1,1:4); n_col = numel(col_names);
for j = 1:n_col
    col_names{j} = sprintf(['<html><center /><font size=5> %s ' ...
        '</font></html>'],char(col_names(j)));
end
row_names = Fx(2:K+1,1);
for i = 1:K
    row_names{i} = sprintf(['<html><font size=5> %s ' ...
        '</font></html>'],char(row_names{i,1}));
end
% now justify content of each cell
for i = 2:K+1
    j = 1; Fx(i,j) = strcat({'         '},num2str(Fx{i,j}),{' '});
    for j = 2:3
        Fx(i,j) = strcat({'       '},num2str(Fx{i,j}),{' '});
    end
    j = 4; Fx(i,j) = strcat({'     '},num2str(Fx{i,j}),{' '});
end
% approximate height of table
ht = min(0.8,(K - 1) * 0.04 + 0.15);
% open new figure called Table 1: Emulator performance training data set
axpos_fig = [ 0.15 0.05 0.4 ht ]; 
% scale in x and y
axpos_fig([1 3]) = axpos_fig([1 3]) * x_mult;
axpos_fig([2 4]) = axpos_fig([2 4]) * y_mult;
% end scale in x and y
% % fig = figure('units','normalized','name',['Table 1: Emulator ' ...
% %     'performance training data set'],'numbertitle','off',...
% %     'outerposition',axpos_fig);
fig = figure('units','normalized','name',['Table 1: ' ...
    'Emulator performance training data set'],'numbertitle','off',...
    'position',axpos_fig);
% plot table
axpos_table = [ 0.1 0.1 0.8 0.73 ];
% scale in x and y
axpos_table([1 3]) = axpos_table([1 3]) * x_mult;
axpos_table([2 4]) = axpos_table([2 4]) * y_mult;
% end scale in x and y
tbl = uitable(fig,'units','normalized','position',axpos_table,...
    'data',Fx(2:K+1,1:4),'columnname',col_names,'rowname', ...
    row_names,'fontsize',fontsize_table);
set(tbl,'columnWidth', {120})
% set row_header/cell_size
set_row_cells(tbl,col_names,row_names);

% Add title
axpos_fig = [ 0.07 0.85 (1 - 2*0.07) 0.12 ];
% scale in x and y
axpos_fig([1 3]) = axpos_fig([1 3]) * x_mult;
axpos_fig([2 4]) = axpos_fig([2 4]) * y_mult;
% end scale in x and y
if K == 1
    uicontrol(fig,'units','normalized','style', 'text', 'position', ...
        axpos_fig , 'string', ...
        ['Table 1: Properties of HDMR emulator, y = f(x), for ' ...
        'training data set (no bootstrapping)'],...
        'fontsize',fontsize_xylabel);
elseif K > 1
    str = strcat(['Table 1: Properties of HDMR emulator, y = f(x), ' ...
        'for randomized training data set of'],...
        {' '},num2str(K),{' '},'bootstrap trials');
    uicontrol(fig,'units','normalized','style', 'text', 'position', ...
        axpos_fig, 'string', str, 'fontsize',fontsize_xylabel);
end
% Save as a PDF
print(gcf,table1_name(id_pdf),'-dpdf','-fillpage'); 
% exportapp(fig,'test.png');

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%                     Plot SA_sig to table in figure                      %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Now print figure with Tabulated results ( similar to HDMR_results.txt )
n_col = size(SA_sig,2); col_names = SA_sig(1,2:n_col);
for j = 1:n_col-1
    col_names{j} = sprintf(['<html><center /><font size=5> %s ' ...
        '</font></html>'],char(col_names(j)));
end
n_row = size(SA_sig,1); row_names = cell(n_row-1,1);
for i = 2:n_row
    row_names{i-1} = sprintf(['<html><font size=5> %s ' ...
        '</font></html>'],char(SA_sig{i,1}));
end
% now justify content of each cell
for i = 2:n_row
    j = 2;  SA_sig(i,j) = strcat({'        '},num2str(SA_sig{i,j}),{' '});
    for j = 3:12
        SA_sig(i,j) = strcat({'     '},SA_sig(i,j),{' '});
    end
    j = 13;  SA_sig(i,j) = strcat({'       '},num2str(SA_sig{i,j}),{' '});
end
% approximate height of table
ht = min(0.8, (n_row - 1) * 0.038 + 0.15); 

% open new figure called Table 1: Emulator performance training data set
axpos_fig = [ 0.02 0.2 0.85 ht ]; 
% scale in x and y
axpos_fig([1 3]) = axpos_fig([1 3]) * x_mult;
axpos_fig([2 4]) = axpos_fig([2 4]) * y_mult;
% end scale in x and y

% Table 2: Variance-based decomposition and sensitivity coefficients"
% fig = figure('units','normalized','name',['Table 3: Variance-based ' ...
%     'decomposition and sensitivity coefficients'],'numbertitle','off',...
%     'outerposition',axpos_fig);
fig = figure('units','normalized','name',['Table 2: Variance-based ' ...
    'decomposition and sensitivity coefficients'],'numbertitle','off',...
    'position',axpos_fig);
% plot table
axpos_table = [ 0.04 0.1 0.90 0.73 ];
% scale in x and y
axpos_table([1 3]) = axpos_table([1 3]) * x_mult;
axpos_table([2 4]) = axpos_table([2 4]) * y_mult;
% end scale in x and y
tbl = uitable(fig,'units','normalized','position',axpos_table,...
    'data',SA_sig(2:n_row,2:n_col),'columnname',col_names, ...
    'rowname',row_names,'fontsize',fontsize_table);
% % set(tbl,'columnWidth', {100,120,100,120,100,120,100,120,100,135,100, ...
% %     120,135,100})
set(tbl,'columnWidth', {100,120,100,120,100,120,100,120,100, ...
    135,100,120,135})
% set row_header/cell_size
set_row_cells(tbl,col_names,row_names)

% Now print title
% get(gcf, 'Position');
axpos_table = [ 0.07 0.85 (1 - 2*0.07) 0.12 ]; 
% scale in x and y
axpos_table([1 3]) = axpos_table([1 3]) * x_mult;
axpos_table([2 4]) = axpos_table([2 4]) * y_mult;
% end scale in x and y

if K == 1
    uicontrol(fig,'units','normalized','style', 'text', 'position', ...
        axpos_table , 'string', ...
        ['Table 2: HDMR results for significant model components ' ...
        'only using all X and Y data (no bootstrapping)'],...
        'fontsize',fontsize_xylabel);
elseif K > 1
    str = strcat(['Table 2: HDMR results for significant model ' ...
        'components only using'],...
        {' '},num2str(K),{' '},'bootstrap trials');
    uicontrol(fig,'units','normalized','style', 'text', 'position', ...
        axpos_table, 'string', str, 'fontsize',fontsize_xylabel);
end
% Add Table to PDF file
print(gcf,table2_name(id_pdf),'-dpdf','-fillpage'); 

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%               Save Figures & Tables to a single PDF file                %
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% Merge figures & tables
input_files = {file_name,table1_name,table2_name};
% Current directory
HDMR_dir = pwd; 
% Forward slash or backward slash?
if contains(HDMR_dir,'/'), slh = '/'; else, slh = '\'; end
% How many PDFs?
n_pdfs = numel(input_files);
% Cell string
input_pdfs = cell(1,n_pdfs);
% Full names of input files
for i = 1:n_pdfs
    input_pdfs(i) = {[HDMR_dir,slh,input_files{i}]};
end
% Name of output file
output_pdf = [HDMR_dir,slh,file_name]; 
% Benjamin Großmann (2021). Merge PDF-Documents 
% (https://www.mathworks.com/matlabcentral/fileexchange/89127- ...
%     merge-pdf-documents), MATLAB Central File Exchange. August 24, 2021.
memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
merger = org.apache.pdfbox.multipdf.PDFMergerUtility;
cellfun(@(f) merger.addSource(f), input_pdfs)
merger.setDestinationFileName(output_pdf)
merger.mergeDocuments(memSet)

% Delete Tables
delete(table1_name); delete(table2_name); 
 
% Open PDF document
open(file_name);

% Print wait statement to the screen
fprintf(' DONE\n');

end
