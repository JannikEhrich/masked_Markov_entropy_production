function [Sigma_min] = est_EP_min_2_2(A12, P2, P3, c3_ini, c4_ini)
%EST_EP_min_2_2 calculates the estimate of the minimum possible entropy
%production for a network with two observable states and two hidden states 
%given the 2- and 3-jump probabilities and the columns sums of the hidden 
%transition probabilities. It requires valid starting guesses for the free 
%parameters c3 and c4
%
% INPUTS:    
%      A12: 2x2 matrix of observed direct jump probabilities
%       P2: 2x2 matrix of observed 2-jump probabilities
%       P3: 2x2 matrix of observed 3-jump probabilities
%   c3_ini: valid starting guess of 1st column sum of hidden matrix
%   c4_ini: valid starting guess of 2nd column sum of hidden matrix
%
% OUTPUTS:  
%   Sigma_min: minimum possible EP
%
% author:   JEhrich
% version:  1.1 (2021-04-12)
% changes:  deleted output of maximum possible EP

global Sigma_ini;
Sigma_ini = est_EP_2_2(A12,P2,P3,c3_ini,c4_ini);

% find minimum
fun = @(c) support_find_min(A12,P2,P3,c(1),c(2));
c_min = fminsearch(fun,[c3_ini,c4_ini]);
Sigma_min = est_EP_2_2(A12,P2,P3,c_min(1),c_min(2));

end


%% support function for finding minimum by treating nan as inf
function Sigma_out = support_find_min(A12,P2,P3,c3,c4)
global Sigma_ini;

Sigma_out = est_EP_2_2(A12,P2,P3,c3,c4);
if isnan(Sigma_out)
    Sigma_out = Sigma_ini;
end

end



