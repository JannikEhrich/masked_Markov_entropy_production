function [Sigma_min,Sigma_max] = est_EP_range_2_2(A12, P2, P3, c3_ini, c4_ini)
%EST_EP_range_2_2 calculates the estimate of the minimum possible and 
%maximum possible entropy production for a network with two observable
%states and two hidden states given the 2- and 3-jump probabilities and the
%columns sums of the hidden transition probabilities. It requires valid
%starting guesses for the free parameters c3 and c4
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
%   Sigma_max: maximum possible EP
%
% author:   JEhrich
% version:  1.0 (2021-03-24)
% changes:  -

% initial EP
Sigma = est_EP_2_2(A12,P2,P3,c3_ini,c4_ini);

global Sigma_min;
Sigma_min = Sigma
global Sigma_max;
Sigma_max = Sigma;

% find minimum
fun = @(c) support_find_min(A12,P2,P3,c(1),c(2));
c_min = fminsearch(fun,[c3_ini,c4_ini]);
Sigma_min = est_EP_2_2(A12,P2,P3,c_min(1),c_min(2));

% find maximum
fun = @(c) support_find_max(A12,P2,P3,c(1),c(2));
c_max = fminsearch(fun,[c3_ini,c4_ini]);
Sigma_max = est_EP_2_2(A12,P2,P3,c_max(1),c_max(2));

end

%% support function for finding minimum by trating nan as inf
function Sigma = support_find_min(A12,P2,P3,c3,c4)
global Sigma_max;

Sigma = est_EP_2_2(A12,P2,P3,c3,c4);
if Sigma > Sigma_max
    Sigma_max = Sigma;
end
if isnan(Sigma)
    Sigma = Sigma_max;
end
end

%% support function for finding maximum by trating nan as -inf
function Sigma = support_find_max(A12,P2,P3,c3,c4)
global Sigma_min;

Sigma = est_EP_2_2(A12,P2,P3,c3,c4);
if Sigma < Sigma_min
    Sigma_min = Sigma;
end
if isnan(Sigma)
    Sigma = Sigma_min;
end
end

