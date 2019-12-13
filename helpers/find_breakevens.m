function [K, k, B, action_order] = find_breakevens(K, k)
% Finds the points at which a linear utility function defined by
% u_d(theta) = K + k * theta is equal for two d (the break-even points)
% Input:
% K: [1 x D] vector of intercepts
% k: [1 x D] vector of gradients
%
% Output:
% K: [1 x D] vector of intercepts ordered such that k_1 < k2 < ...
% k: [1 x D] vector of intercepts ordered such that k_1 < k2 < ...
% B: [1 x D-1] vector of break-even points
% action_order: permutation index to sort decisions such that k_1 < k2 < ...
n_actions = length(k);
[k, action_order] = sort(k);
K = K(action_order);

B = zeros(1, n_actions-1);

for i = 1:(n_actions-1)
    
    B(i) = (K(i) - K(i+1)) / (k(i+1) - k(i));

end

end

