%RANDOM_NETWORKS_INFERENCE_COMPARISON performs and compares the inference 
%of dissipation for random networks with two hidden states
%
% OUTPUTS:  
%       creates eps-figure of comparison
%
% author:   JEhrich
% version:  1.0 (2021-04-14)
% changes:  

clear
close 'all'
clc

%% parameters
% seed
rng(2)

% number of networks
N = 50;
% maximum number of jumps for analytical calculations
n_max = 100;
% accuracy for fitting procedure
accuracy = 1E-4;

%% main loop
% data structures for entropy productions
Sigma = nan(N,1);
Sigma_cg = nan(N,1);
Sigma_DKL = nan(N,1);
Sigma_est = nan(N,1);

% generate random transition matrices sequentially for reproducability
A_mat = nan(4,4,N);
for ii = 1:N
    A_mat(:,:,ii) = gen_random_transition_matrix(ones(4));
end

tic
parfor ii = 1:N
    ii
    A = A_mat(:,:,ii);
    
    %% calculate real entropy production
    p = calc_steady_state(A);
    Sigma(ii) = calc_entropy_production(A,p);
    
    %% calculate coarse-grained EP
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
    p_j = nan(2,2,n_max);
    % sub-matrices
    A12 = A(1:2,1:2);
    B = A(3:4,1:2);
    C = A(1:2,3:4);
    H = A(3:4,3:4);
    % n = 1 are direct jumps
    p_j(:,:,1) = A12;
    % calculate probabilities of jumps with hidden states
    for jj = 2:n_max
        p_j(:,:,jj) = C * H^(jj-2) * B;
    end
    % calculate irreversibility
    Sigma_DKL(ii) = calc_time_series_irr(p_j(:,:,:), p(1:2));
    
    %% estimates from fitting
    % run minization
    Sigma_est(ii) = est_EP_min_2_2(A12,p_j(:,:,2),p_j(:,:,3),...
        A(3,3)+A(4,3),A(3,4)+A(4,4),accuracy);
    
end
toc


%% plotting
% set font size, line width, and marker size
fS = 18;
lW = 2;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

figure();
% plots for legends
loglog(nan,nan,'k','lineWidth',lW,'MarkerSize',mS);
hold on;
loglog(nan,nan,'rx','lineWidth',lW,'MarkerSize',mS);
loglog(nan,nan,'bs','lineWidth',lW,'MarkerSize',mS); 
loglog(nan,nan,'go','lineWidth',lW,'MarkerSize',mS); 
% real plots
loglog(Sigma,Sigma,'k','lineWidth',lW,'MarkerSize',mS);
loglog(Sigma,Sigma_cg,'go','lineWidth',lW,'MarkerSize',mS);
loglog(Sigma,Sigma_DKL,'bs','lineWidth',lW,'MarkerSize',mS); 
loglog(Sigma,Sigma_est,'rx','lineWidth',lW,'MarkerSize',mS); 

axis([min(Sigma),max(Sigma),1E-4,1.5E0]);
xlabel('$\Delta\Sigma$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{est}$',...
    '$\Delta\Sigma_\mathrm{DKL}$','$\Delta\tilde\Sigma$'},...
    'Location','NorthWest');
saveas(gcf, '../doc/random_networks_EP','epsc')









