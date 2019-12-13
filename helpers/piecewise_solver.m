function [n_opt] = piecewise_solver(n_studies, evsi_analytical_ns, c, budget, n_pieces, ub)
% Solves the optimal sample size allocation problem using SOS2 constraints 
% Input:
% n_studies: the number of studies
% evsi_analytical_ns: a function that takes in the sample size and study index and returns the EVSI for the study
% c: a [1 x S] vector of marginal costs for each study (marginal costs are assumed constant)
% budget: the total budget
% n_pieces: the number of piecewise linear pieces to approximate the EVSI function with 
% ub: the upper bound for sample size

% Output:
% n_opt: the optimal sample size allocation
%
% Requires Gurobi (tested with version 8.1.1 and MATLAB 2019b Update 2)
n_grid = linspace(0, ub, n_pieces);
f_values = zeros(n_pieces, n_studies);
for s = 1:n_studies
   for i = 1:n_pieces
    f_values(i, s) = evsi_analytical_ns(n_grid(i), s);
   end
end
n_cell = repmat({n_grid}, 1, n_studies);
one_cell = repmat({ones(1, n_pieces)}, 1, n_studies);

% variable order:
% n_1, ...n n_S, lambda_11, ..., lambda_1n_pieces, ...lambda_Sn_pieces
model.ub = [ones(1, n_studies) * ub ones(1, n_studies*n_pieces) * Inf];
model.obj = [-c reshape(f_values, 1, n_studies*n_pieces)];
model.modelsense = 'Max';
model.A = sparse([c zeros(1, n_studies*n_pieces); % budget
                  zeros(n_studies, n_studies) blkdiag(one_cell{:}); % sum to 1
                  -eye(n_studies), blkdiag(n_cell{:})]); % definition of n
model.rhs = [budget; ones(n_studies, 1); zeros(n_studies, 1)];
model.sense = ['<', repmat('=', 1, n_studies), repmat('=', 1, n_studies)]';

% sos type 2 constraints
for s = 1:n_studies
    model.sos(s).type = 2;
    model.sos(s).index = n_studies+n_pieces*(s-1)+1:n_studies+n_pieces*s;
    model.sos(s).weight = f_values(:, s);
end

result = gurobi(model);

n_opt = result.x(1:n_studies);


end

