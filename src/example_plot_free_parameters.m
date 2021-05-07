%EXAMPLE_PLOT_FREE_PARAMETERS plots the entropy production compatible with
%observed statistics as function of the free parameters for a network with
%2 hidden states
%
% OUTPUTS:  
%       creates eps-figure of heatmap-plot of entropy production rates
%
% author:   JEhrich
% version:  1.2 (2021-05-06)
% changes:  changed to a31 and a42 as free parameters

clear
close 'all'
clc

%% parameters
Dmu = -1;

% number of gridpoints in c3-c4 grid
n_grid = 500;

%% compute observables: 2-jump and 3-jump probabilities
% define transition matrix
A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
     0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
     0.1, 0.1, 0.4, 0.6;
     0.5, 0, 0.2, 0.1];

a_cg_H1 = A(3,1)+A(4,1);
a_cg_H2 = A(3,2)+A(4,2); 
P2 = A(1:2,3:4)*A(3:4,1:2);
P3 = A(1:2,3:4)*A(3:4,3:4)*A(3:4,1:2);


%% calculate real entropy production
% stationary probabilities
p = calc_steady_state(A); 
% calculate real entropy production rate
Sigma = calc_entropy_production(A, p);

%% solution for all combinations
a31_vec = linspace(0,a_cg_H1,n_grid);
a42_vec = linspace(0,a_cg_H2,n_grid);
Sigma_est = nan(n_grid,n_grid);
A_est = nan(4,4,n_grid,n_grid);

tic
for ii = 1:length(a31_vec)
    ii
    a31 = a31_vec(ii);
    for jj = 1:length(a42_vec)
        a42 = a42_vec(jj); 
        % assign B matrix
        a41 = a_cg_H1 - a31;
        a32 = a_cg_H2 - a42;
        B = [a31,a32; a41,a42];
        % check that B is non-singular
        if abs(det(B))>1E-10
            % calculate C matrix
            C = P2*B^(-1);
            % calculate H matrix
            H = C^(-1)*P3*B^(-1);
            % assemble transition matrix
            A_est(:,:,ii,jj) = [A(1:2,1:2), C; B, H];
            % check whether solution is valid
            if ~sum(sum(A_est(:,:,ii,jj)<0))
                q = calc_steady_state(A_est(:,:,ii,jj));
                Sigma_est(ii,jj) = calc_entropy_production(A_est(:,:,ii,jj),q);
            end
        end
        
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
imag = imagesc(a31_vec,a42_vec,Sigma_est');
set(gca,'YDir','normal')
colormap(jet);
set(imag,'AlphaData',~isnan(Sigma_est'));
hold on;
plot(A(3,1),A(4,2),'rx','MarkerSize',20,'LineWidth',3);

cb = colorbar;
set(cb,'TickLabelInterpreter','latex')
xlabel('$a_{31}$','Interpreter','latex');
ylabel('$a_{42}$','Interpreter','latex');
set(gca,'FontSize',fS);
title('$\Delta\Sigma_\mathrm{est}$','Interpreter','latex');
saveas(gcf, '../doc/EP_est_free_parameters','epsc')

%% find minimum EP
min_Sigma = min(Sigma_est(:))
[ind_3_min,ind_4_min] = find(min_Sigma == Sigma_est);

% matrix with lowest EP
A_est(:,:,ind_3_min,ind_4_min)


