function [Sigma_min] = est_EP_min_3_2(pj, c2_ini, c3_ini, accuracy, output)
%EST_EP_min_3_2 calculates the estimate of the minimum possible entropy
%production for a network with 1 observable state and two hidden states 
%given the jump probabilities and the columns sums of the hidden 
%transition probabilities. It requires valid starting guesses for the free 
%parameters c3 and c4
%
% INPUTS:    
%         pj: vector of jump probabilities
%     c2_ini: valid starting guess of 1st column sum of hidden matrix
%     c3_ini: valid starting guess of 2nd column sum of hidden matrix
%   accuracy: set accuracy for estimates
%     output: (optional) decide whether output is needed
%
% OUTPUTS:  
%   Sigma_min: minimum possible EP
%
% author:   JEhrich
% version:  1.0 (2021-04-16)
% changes:  

if nargin < 5
    output = 0;
end

% first grid roughly 100 x 100 with initial c3,c4 included
c2_vec = [linspace(0,c2_ini,ceil(c2_ini*100)), linspace(c2_ini,1,ceil(c2_ini*100))];
c3_vec = [linspace(0,c3_ini,ceil(c3_ini*100)), linspace(c3_ini,1,ceil(c3_ini*100))];
% shrink grid until minimum and maximum found Sigma roughly agree
Sigma_max = inf;
Sigma_min = 0;
scale = 1;
while (Sigma_max-Sigma_min)/Sigma_min > accuracy
    [Sigma_min, c2_vec, c3_vec, Sigma_max] = calculate_grid(pj,c2_vec,c3_vec,scale);
    % output scale and variations
    if output
        disp(['scale: ' num2str(scale) ', min EP: ' num2str(Sigma_min) ', max EP: ' num2str(Sigma_max)]);
    end
    scale = scale/10;
end


end


%% support function for finding minimum and shrinking grid by 10
function [Sigma_min, c2_vec, c3_vec, Sigma_max] = calculate_grid(pj,c2_vec,c3_vec,scale)
% calculate grid of EP
Sigma_grid = nan(length(c2_vec),length(c3_vec));
for ii = 1:length(c2_vec)
    for jj = 1:length(c3_vec)
        Sigma_grid(ii,jj) = est_EP_3_2(pj,c2_vec(ii),c3_vec(jj));
    end
end
% find minimum
Sigma_min = min(Sigma_grid(:));
% find position
[i2,i3] = find(Sigma_min == Sigma_grid,1);
c2_min = c2_vec(i2);
c3_min = c3_vec(i3);
% sub-grid around minimum
c2_vec = [linspace(max(0,c2_min-scale/10),c2_min,50) linspace(c2_min,min(1,c2_min+scale/10),50)];
c3_vec = [linspace(max(0,c3_min-scale/10),c3_min,50) linspace(c3_min,min(1,c3_min+scale/10),50)];
% find maximum
Sigma_max = max(Sigma_grid(:));

end




