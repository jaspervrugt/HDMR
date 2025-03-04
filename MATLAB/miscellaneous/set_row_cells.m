function set_row_cells(tbl,col_names,row_names)
% This function adjusts the width of the row cells using html/javascript

hs = '<html><font size="+2">';  % html start
he = '</font></html>';          % html end
cnh = cellfun(@(x)[hs x he],col_names,'uni',false);
rnh = cellfun(@(x)[hs x he],row_names,'uni',false);
set(tbl,'ColumnName',cnh,'RowName',rnh) 

%get the row header
jscroll=findjobj(tbl);
rowHeaderViewport = jscroll.getComponent(4);
rowHeader = rowHeaderViewport.getComponent(0);
height = rowHeader.getSize; %#ok<NASGU>
rowHeader.setSize(80,360);

%resize the row header
newWidth = 100;
rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
height = rowHeader.getHeight;
rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
rowHeader.setSize(newWidth,height);