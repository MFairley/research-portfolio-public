function [cost, grad1, grad2] = linear_cost_n(n, c)

cost = sum(c .* n, 2);

if nargout > 1
    grad1 = ones(size(n, 1), size(n, 2)) .* c;
end

if nargout > 2
    grad2 = zeros(size(n, 1), size(n, 2));
end

end

