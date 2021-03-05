%ENTROPY_ESTIMATES_EXAMPLE_NETWORK computes analytical results for
%coarse-grained entropy production and time-series irreversibility for
%example network
%
% OUTPUTS:  
%       creates eps-figures showing example jump probabilities and plotting
%       the estimates and the real entropy production per time step
%
% author:   JEhrich
% version:  1.4 (2021-03-01)
% changes:  changed error erstimationt to bootstrapping method

clear
close 'all'
clc

%% parameters
% vector of Delta mu's for analytic calculations
Dmu_vec = linspace(-1,3,500)';
% vector of Delta mu's for simulations
Dmu_vec_sim = linspace(min(Dmu_vec),max(Dmu_vec),20)';
% length of simulated trajectory
T = 1E5;

% number of re-samples for bootstrapping error estimate
n_bs = 100;

% add Dmu=0 at beginning
Dmu_vec = [0; Dmu_vec];
Dmu_vec_sim = [0; Dmu_vec_sim];

%% calculate analytic results
% maximum number of jumps for analytical calculations
n_max = 100;

% data structure for jump probabilities
p_j = nan(2,2,n_max,length(Dmu_vec));
% data-structures for entropy production estimates
Sigma = nan(size(Dmu_vec));
Sigma_cg = nan(size(Dmu_vec));
Sigma_DKL = nan(size(Dmu_vec));

tic
for ii = 1:length(Dmu_vec)
    Dmu = Dmu_vec(ii);
    % define transition matrix
    A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
         0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
         0.1, 0.1, 0.4, 0.6;
         0.5, 0, 0.2, 0.1];
 
    % stationary probabilities
    p = calc_steady_state(A); 
    % calculate real entropy production rate
    Sigma(ii) = calc_entropy_production(A, p);

    %% mapping to Markov model
    % Markovian transition matrix
    A_cg = [A(1:2,1:2), [(A(1,3)*p(3)+A(1,4)*p(4))/(p(3)+p(4)); ...
        (A(2,3)*p(3)+A(2,4)*p(4))/(p(3)+p(4))];
            A(3,1)+A(4,1), A(3,2)+A(4,2), ...
        ((A(3,3)+A(4,3))*p(3) + (A(3,4)+A(4,4))*p(4))/(p(3)+p(4))];
    % coarse-grained stationary probabilities
    p_cg = calc_steady_state(A_cg); 

    % calculate coarse-grained entropy production rate
    Sigma_cg(ii) = calc_entropy_production(A_cg,p_cg);

    %% time-series irreversibility
    % sub-matrices
    A12 = A(1:2,1:2);
    B = A(3:4,1:2);
    C = A(1:2,3:4);
    H = A(3:4,3:4);

    % n = 1 are direct jumps
    p_j(:,:,1,ii) = A12;
    % calculate probabilities of jumps with hidden states
    for jj = 2:n_max
        p_j(:,:,jj,ii) = C * H^(jj-2) * B;
    end

    % calculate irreversibility
    Sigma_DKL(ii) = calc_time_series_irr(p_j(:,:,:,ii), p(1:2));
end
toc

%% simulation
Sigma_cg_sim = nan(size(Dmu_vec_sim));
Sigma_cg_sim_quant = nan(length(Dmu_vec_sim),2);
% max number of jump probabilities in simulation
n_max_sim = 20;
% jump probabilities
p_j_sim = nan(2,2,n_max_sim,length(Dmu_vec_sim));
p_j_sim_quant = nan(2,2,n_max_sim,length(Dmu_vec_sim),2);
Sigma_DKL_sim = nan(size(Dmu_vec_sim));
Sigma_DKL_sim_quant = nan(length(Dmu_vec_sim),2);

tic
for ii = 1:length(Dmu_vec_sim)
    round(ii/length(Dmu_vec_sim)*100)
    Dmu = Dmu_vec_sim(ii);
    % define transition matrix
    A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
         0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
         0.1, 0.1, 0.4, 0.6;
         0.5, 0, 0.2, 0.1];
    
    % simulate trajectory
    [x_traj, ~] = sim_masked_traj(A,T);
    
    % estimate transition matrix and stationary probabilities
    [A_cg_sim, p_sim] = est_trans_probs(x_traj);
    
    % calculate coarse-grained entropy prodcution
    Sigma_cg_sim(ii) = calc_entropy_production(A_cg_sim,p_sim);
    
    % estimate jump probabilities
    p_j_sim(:,:,:,ii) = est_jump_probs(x_traj,n_max_sim);
    
    % calculate time-series irreversibility
    Sigma_DKL_sim(ii) = calc_time_series_irr(p_j_sim(:,:,:,ii), p_sim(1:2));
    
    % re-sample to calculate bootstrap errors
    Sigma_cg_bs = nan(n_bs,1);
    p_j_bs = nan(2,2,n_max_sim,n_bs);
    Sigma_DKL_bs = nan(n_bs,1);
    parfor jj = 1:n_bs
        [x_traj_bs, ~] = sim_masked_traj(A,T);
        [A_cg_sim_bs, p_sim_bs] = est_trans_probs(x_traj_bs);
        Sigma_cg_bs(jj) = calc_entropy_production(A_cg_sim_bs,p_sim_bs);
        p_j_bs(:,:,:,jj) = est_jump_probs(x_traj_bs,n_max_sim);
        Sigma_DKL_bs(jj) = calc_time_series_irr(p_j_bs(:,:,:,jj),p_sim_bs(1:2));
    end
    % .05 quantiles
    Sigma_cg_sim_quant(ii,:) = quantile(Sigma_cg_bs,[0.05,0.95]);
    for jj = 1:2
        for kk = 1:2
            for nn = 1:n_max_sim-1
                p_j_sim_quant(jj,kk,nn,ii,:) = quantile(reshape(p_j_bs(jj,kk,nn,:),n_bs,1),[0.05,0.95]);
            end
        end
    end
    Sigma_DKL_sim_quant(ii,:) = quantile(Sigma_DKL_bs,[0.05,0.95]);
end
toc

%% plot example jump probabilities for Dmu = 0;
% set font size, line width, and marker size
fS = 22;
lW = 2;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% assign errors from quantiles
p_j_sim_err_neg = p_j_sim(:,:,:,1)-p_j_sim_quant(:,:,:,1,1);
p_j_sim_err_pos = p_j_sim_quant(:,:,:,1,2)-p_j_sim(:,:,:,1);

% remove errorbars that are smaller than symbol size
p_j_sim_err_neg(p_j_sim_err_neg./p_j_sim(:,:,:,1) < 0.15) = nan;
p_j_sim_err_pos(p_j_sim_err_pos./p_j_sim(:,:,:,1) < 0.15) = nan;

figure();
% plots for legend
semilogy(nan,nan,'-r^','lineWidth',lW);
hold on;
semilogy(nan,nan,'-go','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-bs','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-kx','lineWidth',lW,'MarkerSize',mS);
% actual plots, uncomment for errorbars
errorbar(0:n_max_sim-1,reshape(p_j_sim(1,2,:,1),n_max_sim,1),reshape(p_j_sim_err_neg(1,2,:),n_max_sim,1),reshape(p_j_sim_err_pos(1,2,:),n_max_sim,1),'bs','lineWidth',lW-0.5,'MarkerSize',mS);
errorbar(0:n_max_sim-1,reshape(p_j_sim(2,1,:,1),n_max_sim,1),reshape(p_j_sim_err_neg(2,1,:),n_max_sim,1),reshape(p_j_sim_err_pos(2,1,:),n_max_sim,1),'go','lineWidth',lW-0.5,'MarkerSize',mS);
errorbar(0:n_max_sim-1,reshape(p_j_sim(1,1,:,1),n_max_sim,1),reshape(p_j_sim_err_neg(1,1,:),n_max_sim,1),reshape(p_j_sim_err_pos(1,1,:),n_max_sim,1),'r^','lineWidth',lW-0.5,'MarkerSize',mS);
errorbar(0:n_max_sim-1,reshape(p_j_sim(2,2,:,1),n_max_sim,1),reshape(p_j_sim_err_neg(2,2,:),n_max_sim,1),reshape(p_j_sim_err_pos(2,2,:),n_max_sim,1),'kx','lineWidth',lW-0.5,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j_sim(1,2,:,1),n_max_sim,1),'bs','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j(1,2,1:n_max_sim,1),n_max_sim,1),'b-','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j_sim(2,1,:,1),n_max_sim,1),'go','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j(2,1,1:n_max_sim,1),n_max_sim,1),'g-','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j_sim(1,1,:,1),n_max_sim,1),'r^','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j(1,1,1:n_max_sim,1),n_max_sim,1),'r-','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j_sim(2,2,:,1),n_max_sim,1),'kx','lineWidth',lW,'MarkerSize',mS);
semilogy(0:n_max_sim-1,reshape(p_j(2,2,1:n_max_sim,1),n_max_sim,1),'k-','lineWidth',lW,'MarkerSize',mS);

xlabel('$n$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$p_{11}(n)$', '$p_{21}(n)$', '$p_{12}(n)$', '$p_{22}(n)$'},...
    'Location','NorthEast');
axis([0,8,2E-4,0.8]);
% save figure
saveas(gcf, '../../doc/example_jump_probs','epsc')


%% plot entropy productions
% assign errors from quantiles
Sigma_cg_sim_err_neg = Sigma_cg_sim(2:end)-Sigma_cg_sim_quant(2:end,1);
Sigma_cg_sim_err_pos = Sigma_cg_sim_quant(2:end,2)-Sigma_cg_sim(2:end);
Sigma_DKL_sim_err_neg = Sigma_DKL_sim(2:end)-Sigma_DKL_sim_quant(2:end,1);
Sigma_DKL_sim_err_pos = Sigma_DKL_sim_quant(2:end,2)-Sigma_DKL_sim(2:end);

% remove errorbars that are smaller than symbol size
Sigma_cg_sim_err_neg(Sigma_cg_sim_err_neg./Sigma_cg_sim(2:end) < 0.2) = nan;
Sigma_cg_sim_err_pos(Sigma_cg_sim_err_pos./Sigma_cg_sim(2:end) < 0.2) = nan;
Sigma_DKL_sim_err_neg(Sigma_DKL_sim_err_neg./Sigma_DKL_sim(2:end) < 0.3) = nan;
Sigma_DKL_sim_err_pos(Sigma_DKL_sim_err_pos./Sigma_DKL_sim(2:end) < 0.3) = nan;

figure();
% plots for legend
semilogy(nan,nan,'-k','lineWidth',lW);
hold on;
semilogy(nan,nan,'-gs','lineWidth',lW,'MarkerSize',mS);
    semilogy(nan,nan,'-bo','lineWidth',lW,'MarkerSize',mS);
% actual plots
semilogy(Dmu_vec(2:end),Sigma(2:end),'-k','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_DKL(2:end),'-g','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_cg(2:end),'-b','lineWidth',lW);
semilogy(Dmu_vec_sim(2:end),Sigma_cg_sim(2:end),'bo','lineWidth',lW,'MarkerSize',mS);
errorbar(Dmu_vec_sim(2:end),Sigma_cg_sim(2:end),Sigma_cg_sim_err_neg,Sigma_cg_sim_err_pos,'bo','lineWidth',lW-0.5,'MarkerSize',mS);
semilogy(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),'gs','lineWidth',lW,'MarkerSize',mS);
errorbar(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),Sigma_DKL_sim_err_neg,Sigma_DKL_sim_err_pos,'gs','lineWidth',lW-0.5,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),1E-6,1]);
xlabel('$\Delta\mu$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{DKL}$',...
    '$\Delta\tilde\Sigma$'},'Location','SouthEast');
% save figure
saveas(gcf, '../../doc/example_EP_1','epsc')




