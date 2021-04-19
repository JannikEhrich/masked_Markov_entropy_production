function [Sigma_est,A_est] = est_EP_3_2(pj, c2, c3)
%EST_EP_3_2 calculates the estimate of the entropy production for a network
%with 1 observable state and two hidden states given the jump
%probabilities and the columns sums of the hidden transition probabilities
%
% INPUTS:    
%    pj: vector of jump probabilities
%    c2: column sum of 1st column of hidden matrix
%    c3: column sum of 2nd column of hidden matrix
%
% OUTPUTS:  
%   Sigma_est: estimate of EP (nan if now valid solution exists)
%       A_est: estimated transition probability matrix
%
% author:   JEhrich
% version:  1.0 (2021-04-16)
% changes:  


% estimated transition probabilities
A_est = nan(3,3);

% fill directly observable transition probabilities
A_est(1,1) = pj(1);

% B submatrix
A_est(2,1) = ((c3 - 1)*A_est(1,1) - pj(2) - c3 + 1)/(c2 - c3);
A_est(3,1) = 1 - A_est(1,1) - A_est(2,1);

% C submatrix
A_est(1,2) = 1 - c2;
A_est(1,3) = 1 - c3;

% H submatrix
A_est(3,2) = ((c2^2 - c2)*pj(2)^2 + ((A_est(1,1) - 1)*c2^2 + (-pj(3) - A_est(1,1) + 1)*c2 + pj(3) - pj(4))*pj(2) - pj(3)*(A_est(1,1) - 1)*c2^2 + pj(4)*(A_est(1,1) - 1)*c2 + pj(3)^2 + pj(3)*(A_est(1,1) - 1) - pj(4)*(A_est(1,1) - 1))/((c2 - c3)*(pj(2)^2 + (A_est(1,1) - 1)*pj(2) - pj(3)*(A_est(1,1) - 1)));
A_est(3,3) = (c3*(c2 - 1)*pj(2)^2 + ((-pj(3) + (c2 - 1)*A_est(1,1) - c2 + 1)*c3 + pj(3) - pj(4))*pj(2) - (A_est(1,1) - 1)*(pj(3)*c2 - pj(4))*c3 + pj(3)^2 + pj(3)*(A_est(1,1) - 1) - pj(4)*(A_est(1,1) - 1))/((c2 - c3)*(pj(2)^2 + (A_est(1,1) - 1)*pj(2) - pj(3)*(A_est(1,1) - 1)));
A_est(2,2) = c2 - A_est(3,2);
A_est(2,3) = c3 - A_est(3,3);

% check if transition matrix is allowed
if sum(sum(A_est < 0 | A_est > 1))
    Sigma_est = nan;
else
    % compute steady state
    p_est = calc_steady_state(A_est);
    Sigma_est = calc_entropy_production(A_est, p_est);
end


end

