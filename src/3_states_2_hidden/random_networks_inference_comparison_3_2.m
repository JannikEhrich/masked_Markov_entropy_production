%RANDOM_NETWORKS_INFERENCE_COMPARISON_3_2 performs and compares the inference 
%of dissipation for random networks with 3 states with two hidden states
%
% OUTPUTS:  
%       creates eps-figure of comparison
%
% author:   JEhrich
% version:  1.0 (2021-04-16)
% changes:  

clear
close 'all'
clc
% add path to support functions
addpath('../');

%% parameters
% seed
rng(1)

% number of networks
N = 1E2;
% maximum number of jumps for analytical calculations
n_max = 7;
% accuracy for fitting procedure
accuracy = 1E-4;

%% main loop
% data structures for entropy productions
Sigma = nan(N,1);
Sigma_cg = nan(N,1);
Sigma_DKL = nan(N,1);
Sigma_est = nan(N,1);

% generate random transition matrices sequentially for reproducability
A_mat = nan(3,3,N);
for ii = 1:N
    A_mat(:,:,ii) = gen_random_transition_matrix(ones(3));
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
    A_cg = [A(1,1), (A(1,2)*p(2)+A(1,3)*p(3))/(p(2)+p(3));
            A(2,1)+A(3,1), 1 - (A(1,2)*p(2)+A(1,3)*p(3))/(p(2)+p(3))];
            
    % coarse-grained stationary probabilities
    p_cg = calc_steady_state(A_cg); 
    % calculate coarse-grained entropy production rate
    Sigma_cg(ii) = calc_entropy_production(A_cg,p_cg);
    
    %% time-series irreversibility
    p_j = nan(n_max,1);
    % sub-matrices
    B = A(2:3,1);
    C = A(1,2:3);
    H = A(2:3,2:3);
    % n = 1 are direct jumps
    p_j(1) = A(1,1);
    % calculate probabilities of jumps with hidden states
    for jj = 2:n_max
        p_j(jj) = C * H^(jj-2) * B;
    end
    
    %% estimates from fitting
    % run minimization
    Sigma_est(ii) = est_EP_min_3_2(p_j,A(2,2)+A(3,2),A(2,3)+A(3,3),accuracy);
    
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
plot(nan,nan,'k','lineWidth',lW,'MarkerSize',mS);
hold on;
plot(nan,nan,'rx','lineWidth',lW,'MarkerSize',mS);
plot(nan,nan,'bs','lineWidth',lW,'MarkerSize',mS); 
plot(nan,nan,'go','lineWidth',lW,'MarkerSize',mS); 
% real plots
plot(Sigma,Sigma,'k','lineWidth',lW,'MarkerSize',mS);
plot(Sigma,Sigma_cg,'go','lineWidth',lW,'MarkerSize',mS);
plot(Sigma,Sigma_DKL,'bs','lineWidth',lW,'MarkerSize',mS); 
plot(Sigma,Sigma_est,'rx','lineWidth',lW,'MarkerSize',mS); 

axis([min(Sigma),max(Sigma),1E-4,1.5E0]);
xlabel('$\Delta\Sigma$','Interpreter','latex');
set(gca,'FontSize',fS);
legend({'$\Delta\Sigma$', '$\Delta\Sigma_\mathrm{est}$',...
    '$\Delta\Sigma_\mathrm{DKL}$','$\Delta\tilde\Sigma$'},...
    'Location','NorthWest');
%saveas(gcf, '../doc/random_networks_EP','epsc')









