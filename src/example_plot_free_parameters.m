%EXAMPLE_PLOT_FREE_PARAMETERS plots the entropy production compatible with
%observed statistics as function of the free parameters for a network with
%2 hidden states
%
% OUTPUTS:  
%       creates eps-figure of heatmap-plot of entropy production rates
%
% author:   JEhrich
% version:  1.5 (2021-05-18)
% changes:  Updated Delta mu

clear
close 'all'
clc

%% parameters
Dmu = 0;

% number of gridpoints in c3-c4 grid
n_grid = 500;

%% compute observables: 2-jump and 3-jump probabilities
% define transition matrix
A = [0.4 - 0.1*exp(Dmu/2), 0.2*exp(-Dmu/2), 0.3, 0.3;
     0.1*exp(Dmu/2), 0.9 - 0.2*exp(-Dmu/2), 0.1, 0;
     0.1, 0.1, 0.4, 0.6;
     0.5, 0, 0.2, 0.1];
p_gen = calc_steady_state(A)
J_gen = A.*p_gen' - (A.*p_gen')'


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
Sigma_fit = nan(n_grid,n_grid);
A_fit = nan(4,4,n_grid,n_grid);

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
            A_fit(:,:,ii,jj) = [A(1:2,1:2), C; B, H];
            % check whether solution is valid
            if ~sum(sum(A_fit(:,:,ii,jj)<0))
                q = calc_steady_state(A_fit(:,:,ii,jj));
                Sigma_fit(ii,jj) = calc_entropy_production(A_fit(:,:,ii,jj),q);
            end
        end
        
    end
end
toc

%% example points
a31_1 = 0.58;
a42_1 = 0.095;
[~,ind1] = min(abs(a31_vec-a31_1))
[~,ind2] = min(abs(a42_vec-a42_1))
A_fit(:,:,ind1,ind2)
p1 = calc_steady_state(A_fit(:,:,ind1,ind2))
J1 = A_fit(:,:,ind1,ind2).*p1' - (A_fit(:,:,ind1,ind2).*p1')'

a31_2 = 0.025;
a42_2 = 0.06;
[~,ind1] = min(abs(a31_vec-a31_2))
[~,ind2] = min(abs(a42_vec-a42_2))
A_fit(:,:,ind1,ind2)
p2 = calc_steady_state(A_fit(:,:,ind1,ind2))
J2 = A_fit(:,:,ind1,ind2).*p2' - (A_fit(:,:,ind1,ind2).*p2')'


%% plotting
% set font size, line width, and marker size
fS = 18;
lW = 2;
mS = 11;
% set interpreter to latex
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

figure();
imag = imagesc(a31_vec,a42_vec,Sigma_fit');
set(gca,'YDir','normal')
colormap(jet);
set(imag,'AlphaData',~isnan(Sigma_fit'));
hold on;
plot(A(3,1),A(4,2),'ro','MarkerSize',13,'LineWidth',3);
plot(a31_1,a42_1,'o','Color',[1 1 1]*0.45,'MarkerSize',13,'LineWidth',3);
plot(a31_2,a42_2,'o','Color',[1 1 1]*0.45,'MarkerSize',13,'LineWidth',3);

cb = colorbar;
set(cb,'TickLabelInterpreter','latex')
xlabel('$a_{31}$','Interpreter','latex');
ylabel('$a_{42}$','Interpreter','latex');
set(gca,'FontSize',fS);
title('$\Delta\Sigma_\mathrm{fit}$','Interpreter','latex');
saveas(gcf, '../doc/figure_parameter_sweep/EP_est_free_parameters','epsc')

%% find minimum EP
min_Sigma = min(Sigma_fit(:))
[ind_3_min,ind_4_min] = find(min_Sigma == Sigma_fit);

% matrix with lowest EP
A_fit(:,:,ind_3_min,ind_4_min)





