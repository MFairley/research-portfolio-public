function [evsi, grad1, grad2] = evsi_normal_normal_analytical_x(x, N, mu0, n0, sigma, K, k, B)
% Computes the EVSI, its first and second derivative for one study
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
%
% Outputs:
% evsi: expected value of sample information
% grad1: first derivative of EVSI
% grad2: diagonal of Hessian (off diagonals are 0)
evsi = N * ((sum((K(1:end-1) - K(2:end) + mu0 * (k(1:end-1) - k(2:end))) ...
    .* normcdf(sqrt(n0) / (x * sigma) * (B - mu0)) + ...
    x * sigma / sqrt(n0) * ...
    (k(2:end) - k(1:end-1)) .* normpdf(sqrt(n0) / (x * sigma) * (B - mu0))) + ...
    K(end) + k(end) * mu0 - max(utility_linear(mu0, K, k)))); % include these constants although they do not affect argmax

if nargout > 1
    grad1 = N * sigma / sqrt(n0) * sum((k(2:end) - k(1:end-1)) .* ...
        normpdf(sqrt(n0) / (x * sigma) * (B - mu0)));
end

if nargout > 2
    grad2 = N * sqrt(n0) * (1 / (x.^3 * sigma)) * sum((k(2:end) - k(1:end-1)) .* (B - mu0).^2 .* ...
        normpdf(sqrt(n0) / (x * sigma) * (B - mu0)));
end
    
end
