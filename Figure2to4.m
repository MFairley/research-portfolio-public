% Author: Michael Fairley
% Date: December 11, 2019
% MATLAB 2019b Update 2
% clear all
rng(123) % set seed
id = 'MATLAB:fplot:NotVectorized';
warning('off',id)
n_studies = 5;
% prior
rho = 0.9;
mu0 = [10, 27, 6.5, 22, 30];
n0 = [150, 150, 110, 120, 90];
% data generation
sigma = [17, 16, 16.5, 18, 25];
% utilities for decisions - linear
n_decisions = 3;
K = -1 + 2 * rand(1, n_decisions);
k = -1 + 2 * rand(1, n_decisions);
% manually specificy
K = [350000, 500000, 600000];
k = [10000, 4000, -10000];
c = [10000, 10000, 10000, 10000, 10000]; % cost per sample
C = 100000000; % budget
N = [150000, 50000, 2000, 90000, 80000];
% manually check that there are no dominated actions
[K, k, B, action_order] = find_breakevens(K, k);
u = @(theta) utility_linear(theta, K, k);
% lower and upper bounds for plotting n
n_l_bnd = 0.0;
n_u_bnd = 1000;

% n grids for plotting
h_fine = 0.001;
h_course = 0.05;
n_grid_course = n_l_bnd:h_course:n_u_bnd;
n_grid_fine = n_l_bnd:h_fine:n_u_bnd;

% plot utility curves
n_points = 100;
theta_low = -10;
theta_high = 50;
u_d = zeros(n_points, n_decisions);
theta_range = linspace(theta_low, theta_high, n_points);
for i = 1:n_points
    u_d(i, :) = u(theta_range(i));
end
figure();
plot(theta_range, u_d) 
legend('d1', 'd2', 'd3')
title('Utility Curves for 3-action numerical verification');

% objective function
sense_opt = -1;
enbs_analytical_n = @(n) enbs_subgroups_normal_normal_analytical_n(n, N, mu0, n0, sigma, K, k, B, c, sense_opt);
options = optimoptions('fmincon','SpecifyObjectiveGradient', true);
sense_plot = 1;
enbs_analytical_ns = @(n, s) enbs_normal_normal_analytical_n(n, N(s), mu0(s), n0(s), sigma(s), K, k, B, c(s));
evsi_analytical_ns = @(n, s) evsi_normal_normal_analytical_n(n, N(s), mu0(s), n0(s), sigma(s), K, k, B);
% constraint
% constraint = @(n) deal(sum(c .* n) - C, []);

% optimization
n_opt_uncon = fmincon(enbs_analytical_n, [500, 366, 500, 500, 500], [], [], [], [], zeros(n_studies, 1), [], [], options);

% optimization over budget grid
budget_grid = linspace(0e6, 5e7, 12*10);
budget_plot_points = [120, 41, 39, 25];
% budget_plot_points = [121, 42, 40, 39, 25];
% insert an additional point, which will be the new B3
% budget_grid = [budget_grid(1:38), 1.3669e+07, budget_grid(39:end)];
budget_grid_n_opt = zeros(n_studies, length(budget_grid));
budget_grid_n_opt_pw = zeros(n_studies, length(budget_grid)); % piecewise solution
n_pieces = 10000;
for i = 1:length(budget_grid)
    budget_grid_n_opt_pw(:, i) = piecewise_solver(n_studies, evsi_analytical_ns, c, budget_grid(i), n_pieces, n_u_bnd);
%     budget_grid_n_opt(:, i) = fmincon(enbs_analytical_n, [500, 500, 500, 500, 500], c, budget_grid(i), [], [], zeros(n_studies, 1), [], [], options);
end

% optimization with "whole" studies
n_opt_uncon = piecewise_solver(n_studies, evsi_analytical_ns, c, Inf, n_pieces, n_u_bnd);
uncon_enbs_values = zeros(1, n_studies);
for s = 1:n_studies
    uncon_enbs_values(s) = enbs_analytical_ns(n_opt_uncon(s), s);
end

model.obj = uncon_enbs_values;
model.modelsense = 'Max';
model.A = sparse(c .* n_opt_uncon'); % budget
model.rhs = budget_grid(budget_plot_points(3));
model.sense = '<';
model.vtype = 'B';
result = gurobi(model);
selected_studies_full = result.x;

% table of results
% sample sizes
budget_grid_n_opt_pw(:, budget_plot_points)
n_opt_uncon .* selected_studies_full

% budget used
budget_used = c * budget_grid_n_opt_pw(:, budget_plot_points) 

% budget left over
budget_grid(:, budget_plot_points) - budget_used

% budget left over for whole studies
budget_used_selected = c * (n_opt_uncon .* selected_studies_full)
budget_grid(budget_plot_points(3)) - budget_used_selected

% enbs values for each solution
enbs_values_plot_points = zeros(1, length(budget_plot_points));
for i = 1:length(budget_plot_points)
    enbs_values_plot_points(1, i) = enbs_analytical_n(budget_grid_n_opt_pw(:, budget_plot_points(i))); 
end

% enbs for whole studies
enbs_analytical_n(n_opt_uncon .* selected_studies_full)

% optimal solution if using the same used up budget for whole studies
n_opt_whole_budget = piecewise_solver(n_studies, evsi_analytical_ns, c, budget_used_selected, n_pieces, n_u_bnd);
budget_used_selected - c * n_opt_whole_budget % none remaining
% 
enbs_analytical_n(n_opt_whole_budget)

% shadow prices
% we know that fifth study is always non-zero
p_star = 5;
% this evsi includes Ns multipled
shadow_prices = zeros(4, 1);
for i = 2:length(budget_plot_points)
    [~, grad5] = evsi_analytical_ns(budget_grid_n_opt_pw(p_star, budget_plot_points(i)), p_star);
    shadow_prices(i-1) = grad5 / c(p_star) - 1;
end

% shadow price for extra case
[~, grad5] = evsi_analytical_ns(n_opt_whole_budget(p_star), p_star);
shadow_price_whole = grad5 / c(p_star) - 1;

% get a continuous version of sadhwo price
shadow_prices_continuous = zeros(length(budget_grid), 1);
for i = 1:length(budget_grid)
    % find a non-zero point
    p_star = 1;
    for p = 1:length(N)
        if budget_grid_n_opt_pw(p, i) > 0
            p_star = p;
        end
    end
    % otherwise will default to 1
    [~, grad5] = evsi_analytical_ns(budget_grid_n_opt_pw(p_star, i), p_star);
    shadow_prices_continuous(i) = grad5 / c(p_star) - 1;
end


% plots
% ordering based on sample size at large budget
[~, s_order] = sort(budget_grid_n_opt_pw(:, end));

cm_g = gray(10);
cm = lines(n_studies);

% evsi curves with cost
for j = 1:n_studies
   i = s_order(j);
   figure()
   hold on
   xlabel('Study sample size')
   ylabel('Value')
   axis([0, 1000, 0, 14e6])
   title(sprintf('Study %i', j))
   enbs_analytical_ni = @(n) enbs_analytical_ns(n, i);
   evsi_analytical_ni = @(n) evsi_analytical_ns(n, i);
   f_values = zeros(length(n_grid_course), 1);
   f_values_evsi = zeros(length(n_grid_course), 1);
   for k = 1:length(n_grid_course)
       f_values(k) = enbs_analytical_ni(n_grid_course(k));
       f_values_evsi(k) = evsi_analytical_ni(n_grid_course(k));
   end
    
   plot(n_grid_course, f_values_evsi, '-' ,'Color' , cm(i, :), 'LineWidth', 4)
   % cost
   plot(n_grid_course, c(i) .* n_grid_course, ':' ,'Color' , cm(i, :), 'LineWidth', 4)
   legend('EVSI', 'Cost')

   set(gca,'FontSize',16)
end

% ENBS curves with budget points
for j = 1:n_studies
   i = s_order(j);
   figure()
   hold on
   xlabel('Study sample size')
   ylabel('ENBS')
   axis([0, 800, -1e6, 5.25e6])
   title(sprintf('Study %i', j))
   enbs_analytical_ni = @(n) enbs_analytical_ns(n, i);
    f_values = zeros(length(n_grid_course), 1);
    for k = 1:length(n_grid_course)
        f_values(k) = enbs_analytical_ni(n_grid_course(k));
    end
    plot(n_grid_course, f_values, '-' ,'Color' , cm(i, :), 'LineWidth', 4)

   for p = 1:length(budget_plot_points)
       if p == 1
          plot(budget_grid_n_opt_pw(i, budget_plot_points(p)), enbs_analytical_ni(budget_grid_n_opt_pw(i, budget_plot_points(p))),...,
           '*', 'Color', cm_g(1, :), 'MarkerSize', 15, 'HandleVisibility','off')
       else
            if p == 2
                label = 4;
            elseif p == 3
                label = 3;
            elseif p == 4
                label = 1;
            end
            plot(budget_grid_n_opt_pw(i, budget_plot_points(p)), enbs_analytical_ni(budget_grid_n_opt_pw(i, budget_plot_points(p))),...,
                '.', 'Color', cm_g(1, :), 'MarkerSize', 20, 'HandleVisibility','off')
            text(budget_grid_n_opt_pw(i, budget_plot_points(p)), enbs_analytical_ni(budget_grid_n_opt_pw(i, budget_plot_points(p)))+0.1e6,...,
                sprintf('B%i', label), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','left')
       end
   end
   % also plot the special budget, B2
   plot(n_opt_whole_budget(i), enbs_analytical_ni(n_opt_whole_budget(i)),...,
                '.', 'Color', cm_g(1, :), 'MarkerSize', 20, 'HandleVisibility','off')
   text(n_opt_whole_budget(i), enbs_analytical_ni(n_opt_whole_budget(i))+0.1e6,...,
                sprintf('B%i', 2), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','left')
   
   plot([0, n_u_bnd], [0, 0], '-k', 'HandleVisibility','off')
   set(gca,'FontSize',16)
end

% plot sample sizes vs budget
figure();
hold on
for j = 1:n_studies
    i = s_order(j);
    plot(budget_grid, budget_grid_n_opt_pw(i, :), '-', 'LineWidth', 4, 'Color' , cm(i, :))
end
% legend('Study 1', 'Study 2', 'Study 3', 'Study 4', 'Study 5', 'Location', 'NorthWest')
text(2.5e7, 210, "Study 1", 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(2.5e7, 364, "Study 2", 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(2.5e7, 406, "Study 3", 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(2.5e7, 444, "Study 4", 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(2.5e7, 706, "Study 5", 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
% add vertical lines at budgets
l1 = xline(budget_grid(budget_plot_points(4)),'-','B1', 'LineWidth', 2);
l2 = xline(budget_used_selected,'-','B2', 'LineWidth', 2);
l3 = xline(budget_grid(budget_plot_points(3)),'-','B3', 'LineWidth', 2);
l4 = xline(budget_grid(budget_plot_points(2)),'-','B4', 'LineWidth', 2);
l1.LabelHorizontalAlignment = 'center';
l2.LabelHorizontalAlignment = 'center';
l3.LabelHorizontalAlignment = 'center';
l4.LabelHorizontalAlignment = 'center';
ylabel('Optimal Sample Size')
xlabel('Budget')
axis([0, 2.5e7, 0, 800])
set(gca,'FontSize',16)

% plot of overall value vs budget
enbs_grid_values = zeros(1, size(budget_grid_n_opt_pw, 2));
for i = 1:length(enbs_grid_values)
    enbs_grid_values(i) = -enbs_analytical_n(budget_grid_n_opt_pw(:, i));
end
figure()
plot(budget_grid, enbs_grid_values, '-', 'LineWidth', 4)
xlabel('Budget')
ylabel('Population ENBS')
axis([0, 2.5e7, 0, 1.5e7])
set(gca,'FontSize',16)

% plot of shadow price
figure;
plot(budget_grid, shadow_prices_continuous, '-', 'LineWidth', 4)
xlabel('Budget')
ylabel('Shadow Price')
set(gca,'FontSize',16)

