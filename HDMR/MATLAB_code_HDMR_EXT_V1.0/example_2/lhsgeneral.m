




function correlatedSamples = lhsgeneral(pd,correlation,n)



% Developed by Iman Moazzen
% 2060 Project, IESVic
% University of Victoria, BC, Canada
% Last Update: August 31, 2016
%
%
% This code is developed to generate different realizations of correlated random variables 
% with any probability distribution function using Latin Hypercube Sampling (LHS) 
%
%
% [correlatedSamples, initialSamples]=lhsgeneral(pd,correlation,n)
% input arguments:
% - pd is a cell defined by user including the probability distribution of all variables.
%   Each element in pd is an object representing the probability distribution.
%   To generate the objects, use "makedist" (MATLAB Function, see MATLAB documentation for more information)
% - correlation is a matrix defined by user including the correlation coefficients between random variables
% - n is the number of realizations
% output:
% - correlatedSamples: different realizations of correlated random variables 
%
%
% Example 1:
% pd=cell(1,2);
% pd{1} = makedist('Normal',0,20);
% pd{2} = makedist('Triangular',0,100,150);
% correlation = [1 -0.8;-0.8 1];
% n = 100000;
% correlatedSamples = lhsgeneral(pd,correlation,n);
%
%
% Example 2:
% pd=cell(1,3);
% pd{1} = makedist('Triangular',0,5,10);
% pd{2} = makedist('Normal',-10,1);
% pd{3} = makedist('Uniform',20,40);
% correlation = [1 0.8 0.6;0.8 1 0.4;0.6 0.4 1];
% n = 100000;
% correlatedSamples = lhsgeneral(pd,correlation,n);
%
%
%
% For more information about the author, please visit www.ece.uvic.ca/~imanmoaz
% For more information about the 2060 project, please visit us at https://onlineacademiccommunity.uvic.ca/2060project/ 




l=length(pd);                                                               % number of variables       
RStar=correlation;                                                          % ideal correlation between variables defined by user

x=lhsdesign(n,l,'smooth','off');                                            % generate latin hypercube samples (MATLAB Function, see MATLAB documentation for more information)
independent_samples=zeros(n,l);                                             % preallocation for the matrix

for i=1:l
    prob=x(:,i);
    independent_samples(:,i) = icdf('normal',prob, 0, 1);                   % map latin hypercube samples to values using inverse cumulative distribution functions
end

R=corr(independent_samples);                                                % correlation between the generated values
P=chol(RStar).';                                                            % Cholesky decomposition of the ideal correlation matrix
Q=chol(R).';                                                                % Cholesky decomposition of the actual correlation matrix, 0,
dependent_samples=independent_samples*(P*inv(Q)).';                         % transformation matrix which adds dependency between normal variables 

uniform_dependent_samples=normcdf(dependent_samples);                       % tranforming normal distibution to uniform distribuiton
                                                                            % this transformation preserves the dependency between the variables


for i=1:l
    transformed_samples(:,i)=icdf(pd{i},uniform_dependent_samples(:,i));    % mapping each unifrom variable to the probability distribution defined by user 
end


correlatedSamples=transformed_samples;

fprintf('\n\n')
predefined_correlation=correlation
fprintf('\n\n')
actual_correlation=corr(correlatedSamples)




