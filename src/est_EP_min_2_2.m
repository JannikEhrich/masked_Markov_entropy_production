function [Sigma_min] = est_EP_min_2_2(A12, P2, P3, c3_ini, c4_ini, accuracy, output)
%EST_EP_min_2_2 calculates the estimate of the minimum possible entropy
%production for a network with two observable states and two hidden states 
%given the 2- and 3-jump probabilities and the columns sums of the hidden 
%transition probabilities. It requires valid starting guesses for the free 
%parameters c3 and c4
%
% INPUTS:    
%        A12: 2x2 matrix of observed direct jump probabilities
%         P2: 2x2 matrix of observed 2-jump probabilities
%         P3: 2x2 matrix of observed 3-jump probabilities
%     c3_ini: valid starting guess of 1st column sum of hidden matrix
%     c4_ini: valid starting guess of 2nd column sum of hidden matrix
%   accuracy: set accuracy for estimates
%     output: (optional) decide whether output is needed
%
% OUTPUTS:  
%   Sigma_min: minimum possible EP
%
% author:   JEhrich
% version:  1.2 (2021-04-14)
% changes:  added output functionaility

if nargin < 7
    output = 0;
end

% first grid roughly 100 x 100 with initial c3,c4 included
c3_vec = [linspace(0,c3_ini,ceil(c3_ini*100)), linspace(c3_ini,1,ceil(c3_ini*100))];
c4_vec = [linspace(0,c4_ini,ceil(c4_ini*100)), linspace(c4_ini,1,ceil(c4_ini*100))];
% shrink grid until minimum and maximum found Sigma roughly agree
Sigma_max = inf;
Sigma_min = 0;
scale = 1;
while (Sigma_max-Sigma_min)/Sigma_min > accuracy
    [Sigma_min, c3_vec, c4_vec, Sigma_max] = calculate_grid(A12,P2,P3,c3_vec,c4_vec,scale);
    % output scale and variations
    if output
        disp(['scale: ' num2str(scale) ', min EP: ' num2str(Sigma_min) ', max EP: ' num2str(Sigma_max)]);
    end
    scale = scale/10;
end


end


%% support function for finding minimum and shrinking grid by 10
function [Sigma_min, c3_vec, c4_vec, Sigma_max] = calculate_grid(A12,P2,P3,c3_vec,c4_vec,scale)
% calculate grid of EP
Sigma_grid = nan(length(c3_vec),length(c4_vec));
for ii = 1:length(c3_vec)
    for jj = 1:length(c4_vec)
        Sigma_grid(ii,jj) = est_EP_2_2(A12,P2,P3,c3_vec(ii),c4_vec(jj));
    end
end
% find minimum
Sigma_min = min(Sigma_grid(:));
% find position
[i3,i4] = find(Sigma_min == Sigma_grid,1);
c3_min = c3_vec(i3);
c4_min = c4_vec(i4);
% sub-grid around minimum
c3_vec = [linspace(max(0,c3_min-scale/10),c3_min,50) linspace(c3_min,min(1,c3_min+scale/10),50)];
c4_vec = [linspace(max(0,c4_min-scale/10),c4_min,50) linspace(c4_min,min(1,c4_min+scale/10),50)];
% find maximum
Sigma_max = max(Sigma_grid(:));

end




