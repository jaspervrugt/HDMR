function Y_sim = sacmodel ( x )
% SAC-SMA model: old implementation using C-code with explicit Euler step

persistent data tout maxT

if isempty ( data )
    % Load the data
    load bound.txt;
    % Now define length run period
    maxT = 3717; tout = [ 1 : maxT ];
    % Now assign rainfall data (daily total in mm); mm/day
    data.P = bound(1:maxT,6:9);
    % Now assign PET data (daily total in mm); mm/day
    data.Ep = bound(1:maxT,5);
end

% Define initial states (initially states are assumed almost empty)
y0 = 0.5*[x(1) x(2) x(8) x(9) x(10) (x(1)+x(8))];
% Now run model and compute time evolution states and river streamflow
[tci,~,bfcc] = sacsma(tout,y0,x(1:13),[data.Ep(1:maxT,1) data.P(1:maxT,1:4)]);
% Now routing
[ssf,~] = routingtot(tout,tci,[0 0 0],[x(14) x(14) x(14)]);
% Calculate total 6 hourly
totflow = (bfcc + ssf);
% Now reshape to get total daily flow
flow = reshape(totflow,4,numel(tout)); flow = flow'; Y_sim = sum(flow,2); Y_sim = Y_sim(66:maxT);