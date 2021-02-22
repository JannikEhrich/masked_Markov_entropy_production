%ENTROPY_ESTIMATES_EXAMPLE_NETWORK computes analytical results for
%coarse-grained entropy production and time-series irreversibility for
%example network
%
% OUTPUTS:  
%       creates eps-figure plotting the estimates and the real entropy
%       production
%
% author:   JEhrich
% version:  1.1 (2021-02-22)
% changes:  added estimation of jump probabilities and time-serieies
% irreversibility from simulated data. Added error estimation of jump
% probabilities

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

% add Dmu=0 at beginning
Dmu_vec = [0; Dmu_vec];
Dmu_vec_sim = [0; Dmu_vec_sim];

%% calculate analytic results
% maximum number of jumps
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

%% calculate simulation results
Sigma_cg_sim = nan(size(Dmu_vec_sim));
% max number of jump probabilities in simulation
n_max_sim = 20;
p_j_sim = nan(2,2,n_max_sim,length(Dmu_vec_sim));
Sigma_DKL_sim = nan(size(Dmu_vec_sim));


tic
for ii = 1:length(Dmu_vec_sim)
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
    
    % estimate jump probabilities
    for jj = 0:n_max_sim-1
        for kk = 1:2
            for ll = 1:2
                substr = [ll, repmat(3,1,jj), kk];
                p_j_sim(kk,ll,jj+1,ii) = length(strfind(x_traj',substr));
            end
        end
    end
    % normalize
    p_j_sim(:,:,:,ii) = p_j_sim(:,:,:,ii)./sum(sum(p_j_sim(:,:,:,ii),3),1);
    
    % calculate time-series irreversibility
    Sigma_DKL_sim(ii) = calc_time_series_irr(p_j_sim(:,:,:,ii), p_sim(1:2));
end

% calculate error for jump probabilities
p_j_sim_err = nan(size(p_j_sim));
for ii = 1:n_max_sim
    n_trial = T - ii;
    p_j_sim_err(:,:,ii,:) = sqrt(p_j_sim(:,:,ii,:).*(1-p_j_sim(:,:,ii,:))/n_trial);
end

toc

%% plot example jump probabilities for Dmu = 0;
% set font size, line width, and marker size
fS = 18;
lW = 2;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

figure();
% plots for legend
semilogy(nan,nan,'-r^','lineWidth',lW);
hold on;
semilogy(nan,nan,'-go','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-bs','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-kx','lineWidth',lW,'MarkerSize',mS);
% actual plots, uncomment for errorbars
% errorbar(0:n_max_sim-1,reshape(p_j_sim(1,2,:,1),n_max_sim,1),reshape(p_j_sim_err(1,2,:,1),n_max_sim,1),'bs','lineWidth',lW,'MarkerSize',mS);
% errorbar(0:n_max_sim-1,reshape(p_j_sim(2,1,:,1),n_max_sim,1),reshape(p_j_sim_err(2,1,:,1),n_max_sim,1),'go','lineWidth',lW,'MarkerSize',mS);
% errorbar(0:n_max_sim-1,reshape(p_j_sim(1,1,:,1),n_max_sim,1),reshape(p_j_sim_err(1,1,:,1),n_max_sim,1),'r^','lineWidth',lW,'MarkerSize',mS);
% errorbar(0:n_max_sim-1,reshape(p_j_sim(2,2,:,1),n_max_sim,1),reshape(p_j_sim_err(2,2,:,1),n_max_sim,1),'kx','lineWidth',lW,'MarkerSize',mS);
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


%% plot entropy productions

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
semilogy(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),'gs','lineWidth',lW,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),1E-6,1]);
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{DKL}$',...
    '$\Delta\tilde\Sigma$'},'Location','SouthEast');



