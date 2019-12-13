function [enbs, grad1, grad2] = enbs_subgroups_normal_normal_analytical_x(x, N, mu0, n0, sigma, K, k, B, c)
% Computes the ENBS, its first and second derivative for each
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
% c: [1 x S] vector of marginal costs
%
% Outputs:
% enbs: expected net benefit of sampling
% grad1: first derivative of ENBS
% grad2: diagonal of Hessian (off diagonals are 0)
n_subgroups = length(x);
enbs = 0;
gradtmp1 = zeros(1, n_subgroups);
gradtmp2 = zeros(1, n_subgroups);
for s = 1:n_subgroups
    [tmp, gradtmp1(s), gradtmp2(s)] = enbs_normal_normal_analytical_x(x(s), N(s), mu0(s), n0(s), sigma(s), K, k, B, c(s));
    enbs = enbs + tmp;
end

if nargout > 1
   grad1 = gradtmp1;
end

if nargout > 2
    grad2 = gradtmp2;
end

end

