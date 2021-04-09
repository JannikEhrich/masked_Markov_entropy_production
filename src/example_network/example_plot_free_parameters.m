%EXAMPLE_PLOT_FREE_PARAMETERS plots the entropy production compatible with
%observed statistics as function of the free parameters for a network with
%2 hidden states
%
% OUTPUTS:  
%       creates eps-figure of heatmap-plot of entropy production rates
%
% author:   JEhrich
% version:  1.0 (2021-03-24)
% changes:  

clear
close 'all'
clc

%% parameters
Dmu = -1;

% number of gridpoints in c3-c4 grid
n_grid = 200;
%rng(2)

%% compute observables: 2-jump and 3-jump probabilities
% define transition matrix
A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
     0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
     0.1, 0.1, 0.4, 0.6;
     0.5, 0, 0.2, 0.1];

%A = genRandomTransitionMatrix(ones(4))

a_cg_H1 = A(3,1)+A(4,1);
a_cg_H2 = A(3,2)+A(4,2); 
P2 = A(1:2,3:4)*A(3:4,1:2);
P3 = A(1:2,3:4)*A(3:4,3:4)*A(3:4,1:2);

%% find min and max
%[Sigma_min, Sigma_max] = est_EP_range_2_2(A(1:2,1:2),P2,P3,0.5,0.69)


%% calculate real entropy production
% stationary probabilities
p = calc_steady_state(A); 
% calculate real entropy production rate
Sigma = calc_entropy_production(A, p);

%% solution for all combinations of c3 and c4
%c3_vec = linspace(0,1,n_grid);
%c4_vec = linspace(0,1,n_grid);
c3_vec = linspace(0.38,0.62,n_grid);
c4_vec = linspace(0.682,0.7015,n_grid);
Sigma_est = nan(n_grid,n_grid);
A_est = nan(4,4,n_grid,n_grid);

tic
for ii = 1:length(c3_vec)
    ii
    for jj = 1:length(c4_vec)
        [Sigma_est(ii,jj),A_est(:,:,ii,jj)] = ...
            est_EP_2_2(A(1:2,1:2),P2,P3,c3_vec(ii),c4_vec(jj));        
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
imag = imagesc(c3_vec,c4_vec,Sigma_est');
colormap(jet);
set(imag,'AlphaData',~isnan(Sigma_est'));
hold on;
plot(A(3,3)+A(4,3),A(3,4)+A(4,4),'rx','MarkerSize',20,'LineWidth',3);

cb = colorbar;
set(cb,'TickLabelInterpreter','latex')
set(gca,'YDir','normal');
xlabel('$c_3$','Interpreter','latex');
ylabel('$c_4$','Interpreter','latex');
set(gca,'FontSize',fS);
title('$\Sigma_\mathrm{est}$','Interpreter','latex');
saveas(gcf, '../../doc/EP_est_free_parameters','epsc')

%% find minimum and maximum EP
min_Sigma = min(Sigma_est(:))
[ind_3_min,ind_4_min] = find(min_Sigma == Sigma_est)

max_Sigma = max(Sigma_est(:))
[c3_max,c4_max] = find(max_Sigma == Sigma_est)

% matrix with lowest EP
A_est(:,:,ind_3_min,ind_4_min)

% matrix with highest EP
A_est(:,:,c3_max,c4_max)

