function [evsi, grad1, grad2] = evsi_subgroups_normal_normal_analytical_x(x, N, mu0, n0, sigma, K, k, B)
% Computes the EVSI, its first and second derivative for each
% study/subgroup
% 
% Inputs:
% x: sqrt(n/(n+n0))
% N: [1 x S] vector of population sizes
% mu0: [1 x S] vector of prior means
% n0: [1 x S] vector of prior sample sizes
% sigma: [1 x S] vector of data generating standard deviations
% K: [1 x D] vector of linear utility function intercepts
% k: [1 x D] vector of linear utility function gradients
% B: [1 x D-1] vector of break-even points
%
% Outputs:
% evsi: expected value of sample information
% grad1: first derivative of EVSI
% grad2: diagonal of Hessian (off diagonals are 0)
n_subgroups = length(x);
evsi = 0;
gradtmp1 = zeros(1, n_subgroups);
gradtmp2 = zeros(1, n_subgroups);
for s = 1:n_subgroups
    [tmp, gradtmp1(s), gradtmp2(s)] = evsi_normal_normal_analytical_x(x(s), N(s), mu0(s), n0(s), sigma(s), K, k, B);
    evsi = evsi + tmp;
end

if nargout > 1
   grad1 = gradtmp1;
end

if nargout > 2 
    grad2 = gradtmp2;
end

end

