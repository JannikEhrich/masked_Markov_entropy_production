function [x_traj,z_traj] = sim_masked_traj(A,T)
%SIM_MASKED_TRAJ simulates stationar masked Markov chain specieified by
%transition matrix A for T time steps and returns vector x_traj of visible
%and z_traj of hidden trajectory
%
% INPUTS:    
%   A: transition matrix
%   T: number of time-steps
%
% OUTPUTS:  
%   x_traj: visible trajectory
%   z_traj: hidden trajectory
%
% author:   JEhrich
% version:  1.0 (2021-02-19)
% changes:  -

% compute steady-state
p = calc_steady_state(A);

% data structure for trajectory
z_traj = nan(T+1,1);

% draw from steady-state
z_traj(1) = find(rand < cumsum(p),1);
for t = 1:T
    % draw next state
    z_traj(t+1) = find(rand < cumsum(A(:,z_traj(t))),1);
end

% visible trajectory
x_traj = z_traj;
x_traj(x_traj > 2) = 3;

end

