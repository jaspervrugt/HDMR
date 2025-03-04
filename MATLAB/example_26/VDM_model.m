function [X,Y,SS] = VDM_model(X,n,d,u0,ode_options,dt, ...
    t_end,mV,tbx)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function returns the simulated counts and summary metrics of the   %
% vandermeer (VDM) predator-prey model                                    % 
%                                                                         %
%  SYNOPSIS                                                               %
%   [X,Y,SS] = VDM_model(X,n,d,u0,ode_options,dt,t_end,mV,tbx)            %
%  where                                                                  %
%    X           Nxd matrix of parameter vectors                          %
%    n           lenght of simulation                                     %
%                                                                         %
% Please check the following papers                                       %
%  Massoud, E.C., J. Huisman, E. Benincà, M.C. Dietze, W. Bouten, and     % 
%      J.A. Vrugt (2018), Probing the limits of predictability: data      %
%      assimilation of chaotic dynamics in complex food webs, Ecology     %
%      Letters, 21 (1), pp 93-103, https://doi.org/10.1111/ele.12876      %
%  Benincà, E., K.D. Jöhnk, R. Heerkloss, and J. Huisman (2009), Coupled  %
%      predator–prey oscillations in a chaotic food web, Ecology Letters, %
%      12, pp. 1367-1378, doi:10.1111/j.1461-0248.2009.01391.x            %
%  Vandermeer, J. (1993), Loose coupling of predator-prey cycles:         %
%      entrainment, chaos, and intermittency in the classic MacArthur     %
%      consumer-resource equations, American Naturalist, 141, pp. 687-716 %
%  Vandermeer, J. (2004), Coupled oscillations in food webs: balancing    % 
%      competition and mutualism in simple ecological models. American    %
%      Naturalist, 163, pp. 857-867                                       %
%  Vandermeer, J. (2006), Oscillating populations and biodiversity        % 
%      maintenance, Bioscience, 56, pp. 967-975                           % 
%                                                                         %
%  © Written by Jasper A. Vrugt & Abdullah Sahin                          %
%    University of California Irvine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

N = size(X,1);                                  % # parameter vectors
Y = nan(N,n,4);                                 % initialize model output 
                                                % N vectors, n simulated counts 
                                                % & k = 4 species

% Built-in ode45 function to solve 4 differential equations t = 0 to t_end
switch tbx                                      % N simulations parallel?
    case 0      % sequential
        prev = 0;                               % initialize ct printing
        for i = 1:N
            [~,u] = ode45(@(t,u) ...
                VDM(t,u,X(i,1:d),mV), ...
                (0:dt:t_end),u0,ode_options);
            Y(i,:,:) = reshape(u,1,n,4);        % Each species in own Y
            fprintf(1, repmat('\b',1,prev));    % print progress
            prev = fprintf(['VDM:MODEL ' ...
                'EVALUATIONS: %3.2f%% DONE'],100*(i/N));
        end
    case 1      % parallel evaluation
        poolobj = gcp;                          % Pool already open?
        if isempty(poolobj)                     % Open workers
            poolobj = parpool('local');
        end
        fprintf(['VDM: %d workers are used to evaluate the ' ...
            'vandermeer population model %d times\n'], ...
            poolobj.NumWorkers,N);
        parfor i = 1:N
            [~,u] = ode45(@(t,u) ...
                VDM(t,u,X(i,1:d),mV), ...
                (0:dt:t_end),u0,ode_options);
            Y(i,:,:) = reshape(u,1,n,4);        % Each species in own Y
            fprintf('.')                        % print progress    
        end
        poolobj = gcp('nocreate');              % delete the pool
        delete(poolobj);
end
[~,n,K] = size(Y); SS = nan(N,8,K); ct = 1;     % Summary statistics 
% # peaks tells us about # predator prey cycles and
%   std tells us about dynamics around mean
for k = 1:K                                     % Loop over each species
    for i = 1:N                                 % # peaks of simulation
        [mxL,mxM] = peakfinder(Y(i,1:n,k),0.2);
        if ~isempty(mxM)
            %peak = [mxL mxM];                  % 1*501 values for each i, there are one min values for each k
            %peak_cell{i,:} = num2cell(peak);   % each peak different dimensions?
            dim = size(mxM);                    
            mn_dis = n/(dim(2) + 1);            % mean distance peaks?
            peak_mn = mean(mxM);                % mean peak amplitude?
            peak_mx = max(mxM);                 % maximum peak
            % coherence between two curves??
            % quarter delay ??
            % m = 1: SD
            % m = 3: number of peaks
            % m = 4: average distance between each peaks
            % m = 5: the average amplitude of peaks
            % m = 6: the max of peak
            % for m = 7 - 9:
            % (1) case original:
            % m = 7: beta  = x(1)
            % m = 8: alpha = x(2)
            % (2) case new:
            % m = 7: beta  = X(1)
            % m = 8: alpha = [X(2), X(3)]
            % case original:
            SS(ct,1:8,k) = [std(Y(i,1:n,k),[],2), ...
                mean(Y(i,1:n,k),2),numel(mxL),...
                mn_dis,peak_mn,peak_mx,X(i,1),X(i,2)];
            ct = ct + 1;
        else
            % do nothing
        end
        % case 'new' % metrics + 1
        % X_avg(i) = (X(i,2) + X(i,3))/2
        %     SS(i,1:8,k) = [std(Y(i,1:n,k),[],2),mean(Y(i,1:n,k),2),numel(maxLoc)...
        %     avg_dis,peak_avg,peak_max,X(i,1),X_avg(i)];
    end
end

end

