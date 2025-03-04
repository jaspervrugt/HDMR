function M = permn(V, N)
% PERMN - permutations with repetition
%   Using two input variables V and N, M = PERMN(V,N) returns all
%   permutations of N elements taken from the vector V, with repetitions.
%   V can be any type of array (numbers, cells etc.) and M will be of the
%   same type as V.  
%
%   Examples:
%     M = permn([1 2 3],2) % returns the 9-by-2 matrix:
%              1     1
%              1     2
%              1     3
%              2     1
%              2     2
%              2     3
%              3     1
%              3     2
%              3     3
%
%     M = permn([99 7],4) % returns the 16-by-4 matrix:
%              99     99    99    99
%              99     99    99     7
%              99     99     7    99
%              99     99     7     7
%              ...
%               7      7     7    99
%               7      7     7     7
%
%   NB Matrix sizes increases exponentially at rate (n^N)*N.
%

% Algorithm using for-loops which can be implemented in C or VB
nv = length(V);
M = zeros(nv^N,N); 
for ii = 1:N
    cc = 1;
    for jj = 1:(nv^(ii-1))
        for kk = 1:nv
            for mm = 1:(nv^(N-ii))
                M(cc,ii) = V(kk);
                cc = cc + 1;
            end
        end
    end
end