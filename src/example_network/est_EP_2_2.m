function [Sigma_est,A_est] = est_EP_2_2(A12, P2, P3, c3, c4)
%EST_EP_2_2 calculates the estimate of the entropy production for a network
%with two observable states and two hidden states given the 2- and 3-jump
%probabilities and the columns sums of the hidden transition probabilities
%
% INPUTS:    
%   A12: 2x2 matrix of observed direct jump probabilities
%    P2: 2x2 matrix of observed 2-jump probabilities
%    P3: 2x2 matrix of observed 3-jump probabilities
%   c3: column sum of 1st column of hidden matrix
%   c4: column sum of 2nd column of hidden matrix
%
% OUTPUTS:  
%   Sigma_est: estimate of EP (nan if now valid solution exists)
%       A_est: estimated transition probability matrix
%
% author:   JEhrich
% version:  1.1 (2021-03-26)
% changes:  added output of estimated transition matrix


% estimated transition probabilities
A_est = nan(4,4);
% fill directly observable transition probabilities
A_est(1:2,1:2) = A12;
% calculate necessary coarse-grained EP
a_cg_H1 = 1 - A12(1,1) - A12(2,1);
a_cg_H2 = 1 - A12(1,2) - A12(2,2);
% fill with B-submatrix
A_est(3,1) = (-P2(1, 1) - P2(2, 1) + (-c4 + 1)*a_cg_H1)/(c3 - c4);
A_est(3,2) = (-P2(1, 2) - P2(2, 2) + (-c4 + 1)*a_cg_H2)/(c3 - c4);
A_est(4,1) = (P2(1, 1) + P2(2, 1) + (c3 - 1)*a_cg_H1)/(c3 - c4);
A_est(4,2) = (P2(1, 2) + P2(2, 2) + (c3 - 1)*a_cg_H2)/(c3 - c4);
% fill with C-submatrix
A_est(1,3) = ((P2(2, 2) + (c3 - 1)*a_cg_H2)*P2(1, 1) - (P2(2, 1) + (c3 - 1)*a_cg_H1)*P2(1, 2))/(P2(1, 2)*a_cg_H1 + P2(2, 2)*a_cg_H1 - P2(1, 1)*a_cg_H2 - P2(2, 1)*a_cg_H2);
A_est(1,4) = ((P2(2, 2) + (c4 - 1)*a_cg_H2)*P2(1, 1) - (P2(2, 1) + (c4 - 1)*a_cg_H1)*P2(1, 2))/(P2(1, 2)*a_cg_H1 + P2(2, 2)*a_cg_H1 - P2(1, 1)*a_cg_H2 - P2(2, 1)*a_cg_H2);
A_est(2,3) = ((P2(1, 2) + (c3 - 1)*a_cg_H2)*P2(2, 1) - (P2(1, 1) + (c3 - 1)*a_cg_H1)*P2(2, 2))/(P2(1, 2)*a_cg_H1 + P2(2, 2)*a_cg_H1 - P2(1, 1)*a_cg_H2 - P2(2, 1)*a_cg_H2);
A_est(2,4) = ((P2(1, 2) + (c4 - 1)*a_cg_H2)*P2(2, 1) - (P2(1, 1) + (c4 - 1)*a_cg_H1)*P2(2, 2))/(P2(1, 2)*a_cg_H1 + P2(2, 2)*a_cg_H1 - P2(1, 1)*a_cg_H2 - P2(2, 1)*a_cg_H2);
% fil with H-submatrix
A_est(3,3) = (((P3(2, 2) + P3(1, 2))*P2(2, 2) + P3(2, 2)*a_cg_H2*(c4 - 1))*P2(1, 1)^2 + (((-P3(2, 2) - P3(1, 2))*P2(2, 1) + (-P3(2, 1) - P3(1, 1))*P2(2, 2) - (c4 - 1)*(P3(2, 1)*a_cg_H2 + P3(2, 2)*a_cg_H1))*P2(1, 2) + ((P3(2, 2) + P3(1, 2))*P2(2, 2) - a_cg_H2*(c4 - 1)*(P3(1, 2) - P3(2, 2)))*P2(2, 1) + (-P3(2, 1) - P3(1, 1))*P2(2, 2)^2 + (-a_cg_H2*(c3 - 1)*P3(1, 1) + a_cg_H1*(c3 + c4 - 2)*P3(1, 2) - a_cg_H2*(c3 + c4 - 2)*P3(2, 1) + P3(2, 2)*a_cg_H1*(c3 - 1))*P2(2, 2) - a_cg_H2*(c4 - 1)*(P3(2, 1)*a_cg_H2 - P3(2, 2)*a_cg_H1)*(c3 - 1))*P2(1, 1) + ((P3(2, 1) + P3(1, 1))*P2(2, 1) + P3(2, 1)*a_cg_H1*(c4 - 1))*P2(1, 2)^2 + ((-P3(2, 2) - P3(1, 2))*P2(2, 1)^2 + ((P3(2, 1) + P3(1, 1))*P2(2, 2) + a_cg_H2*(c3 + c4 - 2)*P3(1, 1) - a_cg_H1*(c3 - 1)*P3(1, 2) + a_cg_H2*(c3 - 1)*P3(2, 1) - P3(2, 2)*a_cg_H1*(c3 + c4 - 2))*P2(2, 1) + (c4 - 1)*a_cg_H1*((P3(2, 1) - P3(1, 1))*P2(2, 2) + (P3(2, 1)*a_cg_H2 - P3(2, 2)*a_cg_H1)*(c3 - 1)))*P2(1, 2) + (P2(2, 2)*a_cg_H1 - P2(2, 1)*a_cg_H2)*(P3(1, 2)*P2(2, 1) - P3(1, 1)*P2(2, 2) + (P3(1, 2)*a_cg_H1 - a_cg_H2*P3(1, 1))*(c3 - 1))*(c4 - 1))/((P2(1, 2)*a_cg_H1 + P2(2, 2)*a_cg_H1 - P2(1, 1)*a_cg_H2 - P2(2, 1)*a_cg_H2)*(P2(1, 1)*P2(2, 2) - P2(1, 2)*P2(2, 1))*(c3 - c4));
A_est(3,4) = (((P3(2, 2) + P3(1, 2))*P2(2, 2) + P3(2, 2)*a_cg_H2*(c4 - 1))*P2(1, 1)^2 + (((-P3(2, 2) - P3(1, 2))*P2(2, 1) + (-P3(2, 1) - P3(1, 1))*P2(2, 2) - (c4 - 1)*(P3(2, 1)*a_cg_H2 + P3(2, 2)*a_cg_H1))*P2(1, 2) + ((P3(2, 2) + P3(1, 2))*P2(2, 2) - a_cg_H2*(c4 - 1)*(P3(1, 2) - P3(2, 2)))*P2(2, 1) + (-P3(2, 1) - P3(1, 1))*P2(2, 2)^2 + 2*(c4 - 1)*(P3(1, 2)*a_cg_H1 - P3(2, 1)*a_cg_H2 + P3(2, 2)*a_cg_H1/2 - a_cg_H2*P3(1, 1)/2)*P2(2, 2) - a_cg_H2*(c4 - 1)^2*(P3(2, 1)*a_cg_H2 - P3(2, 2)*a_cg_H1))*P2(1, 1) + ((P3(2, 1) + P3(1, 1))*P2(2, 1) + P3(2, 1)*a_cg_H1*(c4 - 1))*P2(1, 2)^2 + ((-P3(2, 2) - P3(1, 2))*P2(2, 1)^2 + ((P3(2, 1) + P3(1, 1))*P2(2, 2) - (c4 - 1)*(P3(1, 2)*a_cg_H1 - P3(2, 1)*a_cg_H2 + 2*P3(2, 2)*a_cg_H1 - 2*a_cg_H2*P3(1, 1)))*P2(2, 1) + (c4 - 1)*a_cg_H1*((P3(2, 1) - P3(1, 1))*P2(2, 2) + (c4 - 1)*(P3(2, 1)*a_cg_H2 - P3(2, 2)*a_cg_H1)))*P2(1, 2) + (P2(2, 2)*a_cg_H1 - P2(2, 1)*a_cg_H2)*(c4 - 1)*(P3(1, 2)*P2(2, 1) - P3(1, 1)*P2(2, 2) + (c4 - 1)*(P3(1, 2)*a_cg_H1 - a_cg_H2*P3(1, 1))))/((P2(1, 2)*a_cg_H1 + P2(2, 2)*a_cg_H1 - P2(1, 1)*a_cg_H2 - P2(2, 1)*a_cg_H2)*(P2(1, 1)*P2(2, 2) - P2(1, 2)*P2(2, 1))*(c3 - c4));
A_est(4,3) = c3 - A_est(3,3);
A_est(4,4) = c4 - A_est(3,4);

% check if transition matrix is allowed
if sum(sum(A_est < 0 | A_est > 1))
    Sigma_est = nan;
else
    % compute steady state
    p_est = calc_steady_state(A_est);
    Sigma_est = calc_entropy_production(A_est, p_est);
end


end

