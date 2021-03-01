function [p_j] = est_jump_probs(x_traj, n_max)
%EST_TJUMP_PROBS estimates jump probabilities for a trajectory x_traj of a 
%stationary Markov chain
%
% INPUTS:    
%   x_traj: stationary trajectory of Markov chain
%   n_max: maximum number of time steps per jump
%
% OUTPUTS:  
%   p_j: matrix of jump psoribabilites
%
% author:   JEhrich
% version:  1.0 (2021-03-01)
% changes:  -

% estimate jump probabilities
p_j = nan(2,2,n_max);
for ii = 0:n_max-1
    for jj = 1:2
        for kk = 1:2
            substr = [kk, repmat(3,1,ii), jj];
            p_j(jj,kk,ii+1) = length(strfind(x_traj',substr));
        end
    end
end

% normalize
p_j(:,:,:) = p_j(:,:,:)./sum(sum(p_j(:,:,:),3),1);

end

