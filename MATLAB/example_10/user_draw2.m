function X = user_draw2(mu1,mu2,C1,C2,w1,N)
% Draw a Nxd matrix X using mean µ and covariance matrix Σ 
% Written by Jasper A. Vrugt
% University of California Irvine

% Draw N parameter vectors, X = N x d matrix
switch exist('gmdistribution.m','file')
    case {2,5} % built-in functions
        S1 = nan(2,2,2); S1(:,:,1) = eye(2); S2 = S1; 
        S1(:,:,2) = C1; S2(:,:,2) = C2; 
        gm13 = gmdistribution([mu1 mu2],S1,[w1 , 1 - w1]); 
        gm24 = gmdistribution([mu1 mu2],S2,[w1 , 1 - w1]);
        X = [random(gm13,N) random(gm24,N)]; X = X(:,[1 3 2 4]);
%        X(:,[1 3]) = random(gm13,N); X(:,[2 4]) = random(gm24,N);
    otherwise % built-in & own approach 
        X(:,[1 3]) = mvnrnd(mu1,C1,N); X(:,[2 4]) = mvnrnd(mu2,C2,N); 
        Z = rand(N,1); id = Z < alfa; X(id,[1 3]) = mvnrnd(mu1, ...
            eye(2),sum(id));
        Z = rand(N,1); id = Z < alfa; X(id,[2 4]) = mvnrnd(mu1, ...
            eye(2),sum(id));
end

end

% % % Draw samples
% % switch exist('gmdistribution.m','file')
% %     case {2,5} % function exists
% %         [sig13,sig24] = deal(eye(2)); sig13(:,:,2) = C1; sig24(:,:,2) = C2;
% %         gm13 = gmdistribution([mu1 mu2],sig13,[w1 , 1 - w1]);
% %         gm24 = gmdistribution([mu1 mu2],sig24,[w1 , 1 - w1]);
% %         X(:,[1 3]) = random(gm13,N); X(:,[2 4]) = random(gm24,N);
% %     otherwise % use available built-in functions        
% %         X(:,[1 3]) = mvnrnd(mu,Omega_1,N); X(:,[2 4]) = mvnrnd(mu,Omega_2,N); 
% %         Z = rand(N,1); idx = Z < alfa; X(idx,[1 3]) = mvnrnd(m,eye(2),sum(idx));
% %         Z = rand(N,1); idx = Z < alfa; X(idx,[2 4]) = mvnrnd(m,eye(2),sum(idx));
% % end