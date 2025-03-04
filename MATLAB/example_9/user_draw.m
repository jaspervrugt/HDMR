function X = user_draw(mu1,mu2,C1,C2,w1,N,d)
% Draw a Nxd matrix X using mean µ and covariance matrix Σ 
% Written by Jasper A. Vrugt
% University of California Irvine

% Draw N parameter vectors, X = N x d matrix
switch exist('gmdistribution.m','file')
    case {2,5} % built-in functions
        C = nan(d,d,2); C(:,:,1) = C1; C(:,:,2) = C2;
        gm = gmdistribution([mu1 mu2],C,[w1 , 1 - w1]);
        X = random(gm,N);
    otherwise % built-in & own approach        
        X = mvnrnd(mu2,C2,N); 
        Z = rand(N,1); id = Z < alfa; 
        X(id,1:d) = mvnrnd(mu1,C1,sum(id));
end

end
