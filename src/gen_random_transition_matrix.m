function [ A ] = gen_random_transition_matrix( T )
%GEN_RANDOM_TRANSITION_MATRIX generates a random positive K x K transition 
%matrix A with one-sum columns and a topology defined by matrix the binary
%adjacency matrix T
%
% INPUTS:    
%    T: binary topology (or adjacency) matrix
%
% OUTPUTS:  
%   A: random transition matrix
%
% author:   JEhrich
% version:  1.0 (2020-04-14)
%
% changes:  naming changed to underscore

K = length(T);
A = zeros(K,K);
% for each columns
for ii = 1:K
    % determine number n of outgoing connections
    n = sum(T(:,ii));
    % draw n-1 delimeters of [0,1] interval
    rs  = sort(rand(n-1,1));
    % probability is difference between delimiters
    rVec = nan(n,1);
    rVec(1) = rs(1);
    for jj = 2:n-1
        rVec(jj) = rs(jj)-rs(jj-1);
    end
    rVec(n) = 1-rs(n-1);
    A(T(:,ii)==1,ii) = rVec;
end

end

