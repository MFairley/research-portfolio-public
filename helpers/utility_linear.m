function u = utility_linear(theta, K, k)
% Linear utility function
% Inputs:
% K: [1 x D] vector of intercepts
% k: [1 x D] vector of gradients
%
% Output:
% u: [1 x D] vector of utilities for each decision
u = K + theta * k;
end

