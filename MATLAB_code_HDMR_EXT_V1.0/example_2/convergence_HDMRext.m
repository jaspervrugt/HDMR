clc; clear; close all hidden;
% clear workspace and figures

% How many parameters? 
d = 3; 
% Number of samples to be used
N = 8000;
% Latin hypercube
X = Latin(zeros(1,3),ones(1,3),N); % Note: much faster than built-in lhsdesign function
% Model
Y = sin(2*pi*X(:,1) - pi) + 7*(sin(2*pi*X(:,2) - pi)).^2 + 0.1*(2*pi*X(:,3) - pi).^4 .* sin(2*pi*X(:,1) - pi);

m = 1:10;

for i = 1:10
    % Specify HDMR options
    options = struct('graphics',0,'basis',1,'maxorder',2,'m',m(i),'K',1,'R',N,...
        'alfa',0.99,'method',1,'tolopt',1e-3);
    
    % Now run the HDMR toolbox
    [S,Ss,Fx,Em,XY] = HDMR_EXT(X,Y,options);
    
    Y_e = Em.Y_e;
    
    var_rHDMRext(i) = var(Y_e - Y);
    
    SumS_HDMRext(i) = str2double(S{8,7});
end

subplot(1,2,1)
plot(m,var_rHDMRext)
subplot(1,2,2)
plot(m,SumS_HDMRext)
% save('var_rHDMRext.mat','var_rHDMRext')
% save('SumS_HDMRext.mat','SumS_HDMRext')