function [SA,Fx] = construct_tables(SA,nterms,p0,RMSE,alfa,f_ret,K,C2,C3,n1,n2,n3,n)
% This function returns a Table with a) sensitivity, and b) emulator results

% ----------------------------------------------------------------------- %
%                 FIRST POSTPROCESS SENSITIVITY ESTIMATES                 %
% ----------------------------------------------------------------------- %

% Compute average sensitivity values
S_m = mean(SA.S,2); Sa_m = mean(SA.Sa,2); Sb_m = mean(SA.Sb,2);
% Now calculate sum of each statistic
S_sum = sum(SA.S); Sa_sum = sum(SA.Sa); Sb_sum = sum(SA.Sb);
% Now calculate associated std's
Z = @(p) -sqrt(2) * erfcinv(p*2);
% multiplier
m = Z(1 - alfa/2);
% Compute output statistics for Y (variance)
V_em_div_V_y = bsxfun(@rdivide,SA.V_em,SA.V_y); % Generic call (V_em./V_Y OK in R2016)
V_em_rat = mean(V_em_div_V_y,2); V_em_div_V_Y_sum = sum(V_em_div_V_y);

% Compute standard deviation of bootstrap results
if K > 1
    S_range = m * std(SA.S,0,2);
    Sa_range = m * std(SA.Sa,0,2);
    Sb_range = m * std(SA.Sb,0,2);
    V_em_range = m * std(V_em_div_V_y,0,2);
    S_sum_range = m * std(S_sum);
    Sa_sum_range = m * std(Sa_sum);
    Sb_sum_range = m * std(Sb_sum);
    V_em_div_V_Y_sum_range = m * std(V_em_div_V_Y_sum);
else
    [S_range,Sa_range,Sb_range,V_em_range] = deal(nan(n,1));
end

ST = nan(n1,1); ST_range = nan(n1,1);
% Now compute the total sensitivity of each parameter/coefficient
for r = 1:n1
    [ij] = n1 + find(sum(C2 == r,2) == 1); [ijk] = n1 + n2 + find(sum(C3 == r,2) == 1);
    % use all bootstrap trials to determine total sensitivity + range!
    % ST(i,1) = sum(S_m([i ; ij ; ijk],1)); ST_std(i,1) = nan;
    TS = sum(SA.S([r ; ij ; ijk],1:K),1);
    ST(r,1) = mean(TS); ST_range(r,1) = m * std(TS);
end

% ----------------------------------------------------------------------- %
%               NOW CONSTRUCT TABULATED TABLE WITH RESULTS                %
% ----------------------------------------------------------------------- %

% how many rows of this table
nr = n1 + n2 + n3 + 1;
% initialize row_names of Table and f_i ( = order);
row_names = cell(nr,1); f_ord = nan(nr,1);
% now create row_names + f_i
for r = 1:n1
    f_ord(r) = 1; row_names{r} = strcat('x',num2str(r));
end
for i = 1:n2
    r = i+n1; f_ord(r) = 2;
    row_names{r} = strcat('x',num2str(C2(i,1)),'/x',num2str(C2(i,2)));
end
for i = 1:n3
    r = i+n1+n2; f_ord(r) = 3;
    row_names{r} = strcat('x',num2str(C3(i,1)),'/x',num2str(C3(i,2)),'/x',num2str(C3(i,3)));
end
% add as last row name the sum of the previous rows
row_names{nr} = 'sum';
% now create column names
col_names = {'term','order','S_{a}','std.','S_{b}','std.','S','std.','S_{T}','std.','V(i)/V(Y)','std.',...
    '#select'};
% how many columns?
nc = numel(col_names);

% Now reinitialize SA to become a cell matrix with sensitivity analysis results
SA = cell(nr+1,nc); SA(1,1:nc) = col_names;
% first column stores row_names
SA(2:nr+1,1) = row_names;
% Now fill the columns 2 - 8 of T -> [ f_ord S^a std. S^b std. S std. ]
if K == 1, j = 1:2:5; elseif K > 1, j = 1:6; end
for r = 2:nr
    SA(r,2) = { f_ord(r-1) }; data = [ Sa_m(r-1) Sa_range(r-1) Sb_m(r-1) Sb_range(r-1) S_m(r-1) S_range(r-1) ];
    SA(r,j+2) = arrayfun(@(X) sprintf('%4.3f',X),data(j),'UniformOutput',0);
end
% Now fill the columns 9 - 10 of T -> [ ST std. ]
if K == 1, j = 1; elseif K > 1, j = 1:2; end
for r = 2:n1+1
    data = [ ST(r-1) ST_range(r-1) ];
    SA(r,j+8) = arrayfun(@(X) sprintf('%4.3f',X),data(j),'UniformOutput',0);
end
% Now fill the columns 11 - 12 of T -> [ V(i)/V_Y std.]
if K == 1, j = 7; elseif K > 1, j = 7:8; end
for r = 2:nr
    data = [ V_em_rat(r-1) V_em_range(r-1) ];
    SA(r,j+4) = arrayfun(@(X) sprintf('%4.3f',X),data(j-6),'UniformOutput',0);
end
% Now fill column 13 of T -> [ FX.select ]
for r = 2:nr
    SA(r,13) = { f_ret(r-1) };
end

% Now fill the last row of T ( = 'sum' )
if K > 1
    in_T = [ sum(Sa_m) , Sa_sum_range , sum(Sb_m) , Sb_sum_range , sum(S_m) , S_sum_range , nan , ...
        nan , sum(V_em_rat) , V_em_div_V_Y_sum_range , nan ];
elseif K == 1
    in_T = [ sum(Sa_m) , nan , sum(Sb_m) , nan , sum(S_m) , nan , nan , nan , sum(V_em_rat) , ...
        nan , nan ];
end
% Find appropriate columns in last row
j = 2 + find(~isnan(in_T)); SA(nr+1,j) = arrayfun(@(X) sprintf('%4.3f',X),in_T(j-2),'UniformOutput',0);

% ----------------------------------------------------------------------- %
%               NOW CONSTRUCT TABULATED RESULTS EMULATORS                 %
% ----------------------------------------------------------------------- %

% initialize row_names of Table and f_i ( = order);
row_names = cell(K,1); 
% now create emulator number
for r = 1:K
    row_names{r} = strcat(num2str(r));
end
% now create column names
col_names = {'emulator','# terms','# coefs.','RMSE'};
% how many columns?
nc = numel(col_names);
% Define return structure SA with sensitivity analysis results
Fx = cell(K+1,nc); Fx(1,1:nc) = col_names;
% first column stores row_names
Fx(2:K+1,1) = row_names;
% Now fill the columns 2 - 4 of B -> [ # terms, # coefs, RMSE ]
for k = 2:K+1
    Fx(k,2:nc-1) = { nterms(1,k-1) p0(1,k-1) };
end    
% Now change number of significant digits of RMSE column
Fx(2:K+1,nc) = arrayfun(@(X) sprintf('%4.3f',X),RMSE(1,1:K),'UniformOutput',0);