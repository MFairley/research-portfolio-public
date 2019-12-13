function [cost, grad1, grad2] = linear_cost_x(x, n0, c)

n_nu = @(nu) (nu .* n0) ./ (1 - nu);
n_x = @(x) n_nu(x.^2);
cost = sum(c .* n_x(x));

if nargout > 1
    grad1 = 2 * c * n0 * x / (1-x^2)^2;
end

if nargout > 2
    grad2 = 2 * c * n0 * (3*x^2 + 1) / (1-x^2)^3;
end

end


