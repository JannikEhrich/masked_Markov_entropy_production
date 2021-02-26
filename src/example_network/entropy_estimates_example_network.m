%ENTROPY_ESTIMATES_EXAMPLE_NETWORK computes analytical results for
%coarse-grained entropy production and time-series irreversibility for
%example network
%
% OUTPUTS:  
%       creates eps-figures showing example jump probabilities and plotting
%       the estimates and the real entropy production per time step
%
% author:   JEhrich
% version:  1.3 (2021-02-25)
% changes:  added errorbars to time-series-irreversibility estimate and
% added putput of eps

clear
close 'all'
clc

%% parameters
% vector of Delta mu's for analytic calculations
Dmu_vec = linspace(-1,3,200)';
% vector of Delta mu's for simulations
Dmu_vec_sim = linspace(min(Dmu_vec),max(Dmu_vec),20)';
% length of simulated trajectory
T = 1E5;

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
 
    % calculate real entropy production rate
    Sigma(ii) = calc_entropy_production(A);

    %% mapping to Markov model
    % estimate coarse-grained probabilities
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
Sigma_cg_sim_err = nan(size(Dmu_vec_sim));
% max number of jump probabilities in simulation
n_max_sim = 20;
% counts of starting states
N_s = nan(1,2,1,length(Dmu_vec_sim));
% jump probabilities
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
    
    % error of coarse-grained probabilities
    p_sim_err = sqrt(p_sim.*(1-p_sim)/(T+1));
    % error of coarse-grained transition probabilities
    A_cg_sim_err = nan(3,3);
    for jj = 1:3
        A_cg_sim_err(:,jj) = sqrt(A_cg_sim(:,jj).*(1-A_cg_sim(:,jj))/p_sim(jj)/(T+1));
    end

    % calculate coarse-grained entropy prodcution
    Sigma_cg_sim(ii) = calc_entropy_production(A_cg_sim);
    
    % error propagatrion for entropy production estimate
    Sigma_cg_sim_err_tmp = 0;
    % add errors from transition probabilities
    for jj = 1:3
        for kk = 1:3
            ker1 = p_sim(kk)*log(A_cg_sim(jj,kk)/A_cg_sim(kk,jj));
            ker1(isnan(ker1) | isinf(ker1)) = 0;
            ker2 = (A_cg_sim(jj,kk)*p_sim(kk) - A_cg_sim(kk,jj)*p_sim(jj))/A_cg_sim(jj,kk);
            Sigma_cg_sim_err_tmp = Sigma_cg_sim_err_tmp + (ker1 + ker2)^2*A_cg_sim_err(jj,kk)^2;
        end
    end
    % add errors from stationary probabilities
    for jj = 1:3
        for kk = 1:3
            ker = A_cg_sim(kk,jj)*log(A_cg_sim(kk,jj)/A_cg_sim(jj,kk));
            ker(isnan(ker) | isinf(ker)) = 0;
            Sigma_cg_sim_err_tmp = Sigma_cg_sim_err_tmp + (ker)^2*p_sim_err(jj)^2;
        end
    end

    Sigma_cg_sim_err(ii) = sqrt(Sigma_cg_sim_err_tmp);
    
    % estimate jump probabilities
    for jj = 0:n_max_sim-1
        for kk = 1:2
            for ll = 1:2
                substr = [ll, repmat(3,1,jj), kk];
                p_j_sim(kk,ll,jj+1,ii) = length(strfind(x_traj',substr));
            end
        end
    end
    % counts of starting states
    N_s(1,1,1,ii) = sum(sum(p_j_sim(:,1,:,ii),3));
    N_s(1,2,1,ii) = sum(sum(p_j_sim(:,2,:,ii),3));
    
    % normalize
    p_j_sim(:,:,:,ii) = p_j_sim(:,:,:,ii)./sum(sum(p_j_sim(:,:,:,ii),3),1);
    
    % calculate time-series irreversibility
    Sigma_DKL_sim(ii) = calc_time_series_irr(p_j_sim(:,:,:,ii), p_sim(1:2));
end
toc

%% calculate error for jump probabilities
p_j_sim_err = nan(size(p_j_sim));
for ii = 1:n_max_sim
    p_j_sim_err(:,:,ii,:) = sqrt(p_j_sim(:,:,ii,:).*(1-p_j_sim(:,:,ii,:))./N_s);
end

%% calculate error for time-series irreversibility
Sigma_DKL_sim_err = zeros(size(Sigma_DKL_sim));
for ii = 1:length(Dmu_vec_sim)
    S_err = 0;
    for jj = 1:n_max_sim
        % add error of p_12(n)
        ker12 = p_sim(2)*log(p_j_sim(1,2,jj,ii)./p_j_sim(2,1,jj,ii)) + (p_j_sim(1,2,jj,ii)*p_sim(2) - p_j_sim(2,1,jj,ii)*p_sim(1))/p_j_sim(1,2,jj,ii);
        if ~(isnan(ker12) || isinf(ker12))
            S_err = S_err + ker12^2*p_j_sim_err(1,2,jj,ii)^2;
        end
        % add error of p_21(n)
        ker21 = -p_sim(1)*log(p_j_sim(1,2,jj,ii)./p_j_sim(2,1,jj,ii)) - (p_j_sim(1,2,jj,ii)*p_sim(2) - p_j_sim(2,1,jj,ii)*p_sim(1))/p_j_sim(2,1,jj,ii);
        if ~(isnan(ker21) || isinf(ker21))
            S_err = S_err +ker21^2*p_j_sim_err(2,1,jj,ii)^2;
        end
        % add error of p1
        ker1 = -p_j_sim(2,1,jj,ii)*log(p_j_sim(1,2,jj,ii)./p_j_sim(2,1,jj,ii));
        if ~(isnan(ker1) || isinf(ker1))
            S_err = S_err + ker1^2*p_sim_err(1)^2;
        end
        % add error of p2
        ker2 = p_j_sim(1,2,jj,ii)*log(p_j_sim(1,2,jj,ii)./p_j_sim(2,1,jj,ii));
        if ~(isnan(ker2) || isinf(ker2))
            S_err = S_err + ker2^2*p_sim_err(2)^2;
        end
    end
    Sigma_DKL_sim_err(ii) = sqrt(S_err);
end
 

%% plot example jump probabilities for Dmu = 0;
% set font size, line width, and marker size
fS = 18;
lW = 2;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% delete errorbars in plot smaller than symbol size
p_j_sim_err_plot = p_j_sim_err;
p_j_sim_err_plot(p_j_sim_err./p_j_sim < 1.2E-1) = nan;

figure();
% plots for legend
semilogy(nan,nan,'-r^','lineWidth',lW);
hold on;
semilogy(nan,nan,'-go','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-bs','lineWidth',lW,'MarkerSize',mS);
semilogy(nan,nan,'-kx','lineWidth',lW,'MarkerSize',mS);
% actual plots, uncomment for errorbars
errorbar(0:n_max_sim-1,reshape(p_j_sim(1,2,:,1),n_max_sim,1),reshape(p_j_sim_err_plot(1,2,:,1),n_max_sim,1),'bs','lineWidth',lW,'MarkerSize',mS);
errorbar(0:n_max_sim-1,reshape(p_j_sim(2,1,:,1),n_max_sim,1),reshape(p_j_sim_err_plot(2,1,:,1),n_max_sim,1),'go','lineWidth',lW,'MarkerSize',mS);
errorbar(0:n_max_sim-1,reshape(p_j_sim(1,1,:,1),n_max_sim,1),reshape(p_j_sim_err_plot(1,1,:,1),n_max_sim,1),'r^','lineWidth',lW,'MarkerSize',mS);
errorbar(0:n_max_sim-1,reshape(p_j_sim(2,2,:,1),n_max_sim,1),reshape(p_j_sim_err_plot(2,2,:,1),n_max_sim,1),'kx','lineWidth',lW,'MarkerSize',mS);
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
% remove errorbars that are smaller than symbol size
Sigma_cg_sim_err(Sigma_cg_sim_err./Sigma_cg_sim < 0.3) = nan;
Sigma_DKL_sim_err(Sigma_DKL_sim_err./Sigma_DKL_sim < 0.3) = nan;
% set negative errorbars such that they reach zero if going sub-zero
Sigma_cg_sim_err_neg = Sigma_cg_sim_err;
Sigma_cg_sim_err_neg(Sigma_cg_sim_err >= Sigma_cg_sim) = Sigma_cg_sim(Sigma_cg_sim_err >= Sigma_cg_sim)-10^(-8);

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
errorbar(Dmu_vec_sim(2:end),Sigma_cg_sim(2:end),Sigma_cg_sim_err_neg(2:end),Sigma_cg_sim_err(2:end),'bo','lineWidth',lW-0.5,'MarkerSize',mS);
errorbar(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),Sigma_DKL_sim_err(2:end),'gs','lineWidth',lW-0.5,'MarkerSize',mS);
semilogy(Dmu_vec_sim(2:end),Sigma_DKL_sim(2:end),'gs','lineWidth',lW,'MarkerSize',mS);
axis([min(Dmu_vec),max(Dmu_vec),1E-6,1]);
xlabel('$\Delta\mu$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{DKL}$',...
    '$\Delta\tilde\Sigma$'},'Location','SouthEast');
% save figure
saveas(gcf, '../../doc/example_EP_1','epsc')




