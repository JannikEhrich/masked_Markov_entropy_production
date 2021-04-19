%EXAMPLE_PLOT_FREE_PARAMETERS_3_2 plots the entropy production compatible 
%with observed statistics as function of the free parameters for a network
%with 1 visible and 2 hidden states
%
% OUTPUTS:  
%       creates eps-figure of heatmap-plot of entropy production rates
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
Dmu = -1;
%rng(2);

% number of gridpoints in c2-c3 grid
n_grid = 2E2;

%% compute observables: 2-jump and 3-jump probabilities
% define transition matrix
A = [0.4, 0.2, 0.2;
     0.1, 0.2, 0.3;
     0.5, 0.6, 0.5];
sum(A)

%A = gen_random_transition_matrix(ones(3))

% calculate jump probabilities
pj = nan(4,1);
pj(1) = A(1,1);
for ii = 2:length(pj)
    pj(ii) = A(1,2:3)*A(2:3,2:3)^(ii-2)*A(2:3,1);
end
    

%% find min and max
%[Sigma_min, Sigma_max] = est_EP_range_2_2(A(1:2,1:2),P2,P3,0.5,0.69)


%% calculate real entropy production
% stationary probabilities
p = calc_steady_state(A); 
% calculate real entropy production rate
Sigma = calc_entropy_production(A, p);

%% solution for all combinations of c3 and c4
c2_vec = linspace(0,1,n_grid);
c3_vec = linspace(0,1,n_grid);
Sigma_est = nan(n_grid,n_grid);
A_est = nan(3,3,n_grid,n_grid);

tic
for ii = 1:length(c2_vec)
    ii
    for jj = 1:length(c3_vec)
        [Sigma_est(ii,jj),A_est(:,:,ii,jj)] = ...
            est_EP_3_2(pj,c2_vec(ii),c3_vec(jj));        
    end
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
imag = imagesc(c2_vec,c3_vec,Sigma_est');
colormap(jet);
set(imag,'AlphaData',~isnan(Sigma_est'));
hold on;
plot(A(2,2)+A(3,2),A(2,3)+A(3,3),'rx','MarkerSize',20,'LineWidth',3);

cb = colorbar;
set(cb,'TickLabelInterpreter','latex')
set(gca,'YDir','normal');
xlabel('$c_2$','Interpreter','latex');
ylabel('$c_3$','Interpreter','latex');
set(gca,'FontSize',fS);
title('$\Delta\Sigma_\mathrm{est}$','Interpreter','latex');
%saveas(gcf, '../doc/EP_est_free_parameters','epsc')

%% find minimum and maximum EP
% real EP
Sigma

min_Sigma = min(Sigma_est(:))
[ind_2_min,ind_3_min] = find(min_Sigma == Sigma_est);

max_Sigma = max(Sigma_est(:))
[ind_2_max,ind_3_max] = find(max_Sigma == Sigma_est);

% matrix with lowest EP
A_est(:,:,ind_2_min,ind_3_min)

% matrix with highest EP
A_est(:,:,ind_2_max,ind_3_max);

