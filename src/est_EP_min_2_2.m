function [Sigma_min] = est_EP_min_2_2(A12, P2, P3, a31_ini, a42_ini, accuracy, output)
%EST_EP_min_2_2 calculates the estimate of the minimum possible entropy
%production for a network with two observable states and two hidden states 
%given the 2- and 3-jump probabilities and the columns sums of the hidden 
%transition probabilities. It requires valid starting guesses for the free 
%parameters a31 and a42
%
% INPUTS:    
%        A12: 2x2 matrix of observed direct jump probabilities
%         P2: 2x2 matrix of observed 2-jump probabilities
%         P3: 2x2 matrix of observed 3-jump probabilities
%    a31_ini: valid starting guess for first free parameter
%    a42_ini: valid starting guess for second free parameter
%   accuracy: set accuracy for estimates
%     output: (optional) decide whether output is needed
%
% OUTPUTS:  
%   Sigma_min: minimum possible EP
%
% author:   JEhrich
% version:  1.3 (2021-05-06)
% changes:  changed parametrization of solution

if nargin < 7
    output = 0;
end

% first grid roughly 100 x 100 with initial c3,c4 included
a31_vec = [linspace(0,a31_ini,ceil(a31_ini*100)), linspace(a31_ini,1,ceil((1-a31_ini)*100))];
a42_vec = [linspace(0,a42_ini,ceil(a42_ini*100)), linspace(a42_ini,1,ceil((1-a42_ini)*100))];
% shrink grid until minimum and maximum found Sigma roughly agree
Sigma_max = inf;
Sigma_min = 0;
scale = 1;
while (Sigma_max-Sigma_min)/Sigma_min > accuracy
    [Sigma_min, a31_vec, a42_vec, Sigma_max] = calculate_grid(A12,P2,P3,a31_vec,a42_vec,scale);
    % output scale and variations
    if output
        disp(['scale: ' num2str(scale) ', min EP: ' num2str(Sigma_min) ', max EP: ' num2str(Sigma_max)]);
    end
    scale = scale/10;
end


end


%% support function for finding minimum and shrinking grid by 10
function [Sigma_min, a31_vec, a42_vec, Sigma_max] = calculate_grid(A12,P2,P3,a31_vec,a42_vec,scale)
% calculate grid of EP
Sigma_grid = nan(length(a31_vec),length(a42_vec));
a_cg_H1 = 1-A12(1,1)-A12(2,1);
a_cg_H2 = 1-A12(1,2)-A12(2,2);
for ii = 1:length(a31_vec)
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
            A_est = [A12, C; B, H];
            % check whether solution is valid
            if ~sum(sum(A_est<0))
                % sanity-check whether column sums are 1
                if max(abs(sum(A_est)-1)) > 1E-1
                    A_est
                    sum(A_est)
                    error('Column sums are not correct');
                end
                q = calc_steady_state(A_est);
                Sigma_grid(ii,jj) = calc_entropy_production(A_est,q);
            end
        end
        
    end
end
% find minimum
Sigma_min = min(Sigma_grid(:));
% find position
[i1,i2] = find(Sigma_min == Sigma_grid,1);
a31_min = a31_vec(i1);
a42_min = a42_vec(i2);
% sub-grid around minimum
a31_vec = [linspace(max(0,a31_min-scale/10),a31_min,50) linspace(a31_min,min(1,a31_min+scale/10),50)];
a42_vec = [linspace(max(0,a42_min-scale/10),a42_min,50) linspace(a42_min,min(1,a42_min+scale/10),50)];
% find maximum
Sigma_max = max(Sigma_grid(:));

end




