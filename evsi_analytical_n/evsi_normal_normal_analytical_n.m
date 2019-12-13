function [evsi, grad1, grad2] = evsi_normal_normal_analytical_n(n, N, mu0, n0, sigma, K, k, B)
% Computes the EVSI, its first and second derivative for one study
% 
% Inputs:
% n: sample size
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
x = sqrt(n / (n0 + n));
[evsi, grad1x, grad2x] = evsi_normal_normal_analytical_x(x, N, mu0, n0, sigma, K, k, B);

if nargout > 1
    grad_x_n = n0/(2 * sqrt(n / (n0 + n)) * (n0 + n)^2);
    grad1 = grad_x_n * grad1x;
end

if nargout > 2
    grad2_x_n = -(n0 * (n0 + 4 * n))/(4 * (n/(n0 + n))^(3/2) * (n0 + n)^4);
    grad2 = grad2x * (grad_x_n)^2 + grad1x * grad2_x_n;
end

end


