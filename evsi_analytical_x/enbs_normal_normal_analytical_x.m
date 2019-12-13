function [enbs, grad1, grad2] = enbs_normal_normal_analytical_x(x, N, mu0, n0, sigma, K, k, B, c)
% Computes the ENBS, its first and second derivative for one study
% 
% Inputs:
% x: sqrt(n/(n+n0))
% N: population size
% mu0: prior mean
% n0: vector of prior sample size
% sigma: data generating standard deviation
% K: [1 x D] vector of linear utility function intercepts
% k: [1 x D] vector of linear utility function gradients
% B: [1 x D-1] vector of break-even points
% c: marginal cost
%
% Outputs:
% enbs: expected net benefit of sampling
% grad1: first derivative of ENBS
% grad2: diagonal of Hessian (off diagonals are 0)
[evsi, evsi_g1, evsi_g2] = evsi_normal_normal_analytical_x(x, N, mu0, n0, sigma, K, k, B);
[cost, cost_g1, cost_g2] = linear_cost_x(x, n0, c);
    
enbs = evsi - cost;

if nargout > 1
    grad1 = evsi_g1 - cost_g1;
end

if nargout > 2
    grad2 = evsi_g2 - cost_g2;
end
    
end

