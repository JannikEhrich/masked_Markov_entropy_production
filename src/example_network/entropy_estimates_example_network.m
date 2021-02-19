%ENTROPY_ESTIMATES_EXAMPLE_NETWORK computes analytical results for
%coarse-grained entropy production and time-series irreversibility for
%example network
%
% OUTPUTS:  
%       creates eps-figure plotting the estimates and the real entropy
%       production
%
% author:   JEhrich
% version:  1.0 (2021-02-19)
% changes:  

clear
close 'all'
clc

%% parameters
% vector of Delta mu's for analytic calculations
Dmu_vec = linspace(-1,3,500)';
% vector of Delta mu's for simulations
Dmu_vec_sim = linspace(min(Dmu_vec),max(Dmu_vec),30);
% length of simulated trajectory
T = 1E5;

%% calculate analytic results
% data-structures for entropy production estimates
Sigma = nan(size(Dmu_vec));
Sigma_cg = nan(size(Dmu_vec));
Sigma_DKL = nan(size(Dmu_vec));

tic
parfor ii = 1:length(Dmu_vec)
    Dmu = Dmu_vec(ii);
    % define transition matrix
    A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
         0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
         0.1, 0.1, 0.4, 0.6;
         0.5, 0, 0.2, 0.1];
 
    % calculate real entropy production rate
    Sigma(ii) = calc_entropy_production(A);

    %% mapping to Markov model
    % calculate steady-state
    p = calc_steady_state(A);

    % Markovian transition matrix
    A_cg = [A(1:2,1:2), [(A(1,3)*p(3)+A(1,4)*p(4))/(p(3)+p(4)); ...
        (A(2,3)*p(3)+A(2,4)*p(4))/(p(3)+p(4))];
            A(3,1)+A(4,1), A(3,2)+A(4,2), ...
        ((A(3,3)+A(4,3))*p(3) + (A(3,4)+A(4,4))*p(4))/(p(3)+p(4))];

    % calculate coarse-grained entropy production rate
    Sigma_cg(ii) = calc_entropy_production(A_cg);

    %% time-series irreversibility
    % sub-matrices
    A12 = A(1:2,1:2);
    B = A(3:4,1:2);
    C = A(1:2,3:4);
    H = A(3:4,3:4);

    % maximum number of jumps
    n_max = 1E3;
    % matrix of jump probabilities
    p_j = nan(2,2,n_max);
    % n = 1 are direct jumps
    p_j(:,:,1) = A12;
    % calculate probabilities of jumps with hidden states
    for jj = 2:n_max
        p_j(:,:,jj) = C * H^(jj-2) * B;
    end

    % calculate irreversibility
    Sigma_DKL(ii) = calc_time_series_irr(p_j, p(1:2));
end
toc

%% calculate simulation results
Sigma_cg_sim = nan(size(Dmu_vec_sim));
tic
parfor ii = 1:length(Dmu_vec_sim)
    Dmu = Dmu_vec_sim(ii);
    % define transition matrix
    A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
         0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
         0.1, 0.1, 0.4, 0.6;
         0.5, 0, 0.2, 0.1];
    % simulate trajectory
    [x_traj, ~] = sim_masked_traj(A,T);
    % coarse-grained stationary probabilities
    p_sim = [sum(x_traj==1);sum(x_traj==2);sum(x_traj==3)]/(T+1);
    % coarse-grained transition rates
    A_cg_sim = nan(3,3);
    % count transitions
    for jj = 1:3
        for kk = 1:3
            A_cg_sim(jj,kk) = sum(x_traj(2:end)==jj & x_traj(1:end-1)==kk);
        end
    end
    % normalize
    A_cg_sim = A_cg_sim./sum(A_cg_sim);
    
    % calculate coarse-grained entropy prodcution
    Sigma_cg_sim(ii) = calc_entropy_production(A_cg_sim);
end
toc


%% plot entropy productions
% set font size, line width, and marker size
fS = 18;
lW = 2;
mS = 10;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

figure();
% plots for legend
semilogy(nan,nan,'-k','lineWidth',lW);
hold on;
semilogy(nan,nan,'-gs','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-bo','lineWidth',lW,'MarkerSize',mS);
% actual plots
semilogy(Dmu_vec,Sigma,'-k','lineWidth',lW);
semilogy(Dmu_vec,Sigma_DKL,'-g','lineWidth',lW);
semilogy(Dmu_vec,Sigma_cg,'-b','lineWidth',lW);
semilogy(Dmu_vec_sim,Sigma_cg_sim,'bo','lineWidth',lW,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),1E-6,1]);
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{DKL}$',...
    '$\Delta\tilde\Sigma$'},'Location','SouthEast');



