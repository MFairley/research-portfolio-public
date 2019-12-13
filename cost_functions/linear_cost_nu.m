function [cost, grad1, grad2] = linear_cost_nu(nu, n0, c)

n_nu = @(nu) (nu * n0) ./ (1 - nu);
cost = sum(c .* n_nu(nu));

if nargout > 1
    grad1 = c * n0 * (1-nu)^(-2);
end

if nargout > 2
    grad2 = c * n0 * 2/(1 - nu)^3;
end

end