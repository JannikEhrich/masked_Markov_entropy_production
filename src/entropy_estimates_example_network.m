%ENTROPY_ESTIMATES_EXAMPLE_NETWORK computes analytical results for
%coarse-grained entropy production and time-series irreversibility for
%example network
%
% OUTPUTS:  
%       creates eps-figures showing example jump probabilities and plotting
%       the estimates and the real entropy production per time step
%
% author:   JEhrich
% version:  1.9 (2021-07-02)
% changes:  added superscript "min" to "fit"

clear
close 'all'
clc

%% parameters
% vector of Delta mu's for analytic calculations
Dmu_vec = linspace(-1,2.7,100)';
% vector of Delta mu's for simulations
Dmu_vec_sim = linspace(min(Dmu_vec),max(Dmu_vec),20)';
% length of simulated trajectory
T = 1E6;
% accuracy for parameter sweep entropy estimates
accuracy = 1E-4;

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
Sigma_fit = nan(size(Dmu_vec));

tic
for ii = 1:length(Dmu_vec)
    Dmu = Dmu_vec(ii)
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
    
    %% estimates from fitting
    % run minization
    Sigma_fit(ii) = est_EP_min_2_2(A12,p_j(:,:,2,ii),p_j(:,:,3,ii),A(3,1),A(4,2),accuracy);
      
end
toc

%% simulation
Sigma_cg_sim = nan(size(Dmu_vec_sim));
Sigma_cg_sim_quant = nan(length(Dmu_vec_sim),2);
% max number of jump probabilities in simulation
n_max_sim = 20;
% jump probabilities
p_j_sim = nan(2,2,n_max_sim,length(Dmu_vec_sim));
Sigma_DKL_sim = nan(size(Dmu_vec_sim));
Sigma_fit_sim = nan(size(Dmu_vec_sim));

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

    % estimates from fitting
    % run minization
    Sigma_fit_sim(ii) = est_EP_min_2_2(A_cg_sim(1:2,1:2),p_j_sim(:,:,2,ii),p_j_sim(:,:,3,ii),A(3,1),A(4,2),accuracy);

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
% actual plots
semilogy(1:n_max_sim,reshape(p_j_sim(1,2,:,1),n_max_sim,1),'bs','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j(1,2,1:n_max_sim,1),n_max_sim,1),'b-','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j_sim(2,1,:,1),n_max_sim,1),'go','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j(2,1,1:n_max_sim,1),n_max_sim,1),'g-','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j_sim(1,1,:,1),n_max_sim,1),'r^','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j(1,1,1:n_max_sim,1),n_max_sim,1),'r-','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j_sim(2,2,:,1),n_max_sim,1),'kx','lineWidth',lW,'MarkerSize',mS);
semilogy(1:n_max_sim,reshape(p_j(2,2,1:n_max_sim,1),n_max_sim,1),'k-','lineWidth',lW,'MarkerSize',mS);

xlabel('$n$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$p_{11}(n)$', '$p_{21}(n)$', '$p_{12}(n)$', '$p_{22}(n)$'},...
    'Location','NorthEast');
axis([1,9,2E-4,0.8]);
% save figure
saveas(gcf, '../doc/example_jump_probs','epsc')


%% plot entropy productions only coarse-grained estimate
figure();
% plots for legend
semilogy(nan,nan,'-k','lineWidth',lW);
hold on;
semilogy(nan,nan,'-go','lineWidth',lW,'MarkerSize',mS);
% actual plots
semilogy(Dmu_vec(2:end),Sigma(2:end),'-k','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_cg(2:end),'-g','lineWidth',lW);
semilogy(Dmu_vec_sim(2:end),Sigma_cg_sim(2:end),'go','lineWidth',lW,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),3E-6,0.6]);
xlabel('$\Delta\mu$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$','$\Delta\tilde\Sigma$'},'Location','SouthEast');
% save figure
saveas(gcf, '../doc/example_EP_0','epsc')

%% plot entropy productions without Sigma_fit
figure();
% plots for legend
semilogy(nan,nan,'-k','lineWidth',lW);
hold on;
semilogy(nan,nan,'-bs','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-go','lineWidth',lW,'MarkerSize',mS);
% actual plots
semilogy(Dmu_vec(2:end),Sigma(2:end),'-k','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_DKL(2:end),'-b','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_cg(2:end),'-g','lineWidth',lW);
semilogy(Dmu_vec_sim(2:end),Sigma_cg_sim(2:end),'go','lineWidth',lW,'MarkerSize',mS);
semilogy(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),'bs','lineWidth',lW,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),3E-6,0.6]);
xlabel('$\Delta\mu$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{DKL}$',...
    '$\Delta\tilde\Sigma$'},'Location','SouthEast');
% save figure
saveas(gcf, '../doc/example_EP_1','epsc')

%% plot entropy productions with Sigma_fit
figure();
% plots for legend
semilogy(nan,nan,'-k','lineWidth',lW);
hold on;
semilogy(nan,nan,'-rx','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-bs','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-go','lineWidth',lW,'MarkerSize',mS);
% actual plots
semilogy(Dmu_vec(2:end),Sigma(2:end),'-k','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_fit(2:end),'-r','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_DKL(2:end),'-b','lineWidth',lW);
semilogy(Dmu_vec(2:end),Sigma_cg(2:end),'-g','lineWidth',lW);
semilogy(Dmu_vec_sim(2:end),Sigma_cg_sim(2:end),'go','lineWidth',lW,'MarkerSize',mS);
semilogy(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),'bs','lineWidth',lW,'MarkerSize',mS);
semilogy(Dmu_vec_sim(2:end),Sigma_fit_sim(2:end),'rx','lineWidth',lW,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),3E-6,0.6]);
xlabel('$\Delta\mu$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{fit}^\mathrm{min}$',...
    '$\Delta\Sigma_\mathrm{DKL}$','$\Delta\tilde\Sigma$'},...
    'Location','SouthEast');
% save figure
saveas(gcf, '../doc/example_EP_2','epsc')



