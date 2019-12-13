% Author: Michael Fairley
% Date: December 11, 2019
% MATLAB 2019b Update 2
rng(123) % set seed
id = 'MATLAB:fplot:NotVectorized';
warning('off',id)
% prior
mu0 = 10;
n0 = 150;
% data generation
sigma = 17;
% utilities for decisions - linear
% manually specificy
n_decisions1 = 3;
n_decisions2 = 2;
K1 = [350000, 500000, 600000];
k1 = [10000, 4000, -10000];
K2 = [450000, 500000];
k2 = [10000, 4000];
c = 10000; % cost per sample
C = 100000000; % budget
N = [100000, 15000];
% manually check that there are no dominated actions
[K1, k1, B1, action_order1] = find_breakevens(K1, k1);
[K2, k2, B2, action_order2] = find_breakevens(K2, k2);
u = @(theta, K, k) utility_linear(theta, K, k);
% lower and upper bounds for plotting n
n_l_bnd = 0.0;
n_u_bnd = 1000;

% n grids for plotting
h_fine = 0.001;
h_course = 0.05;
n_grid_course = n_l_bnd:h_course:n_u_bnd;
n_points_n = length(n_grid_course);


% plot utility curves
n_points = 100;
theta_low = -10;
theta_high = 50;
u_d1 = zeros(n_points, n_decisions1);
u_d2 = zeros(n_points, n_decisions2);
theta_range = linspace(theta_low, theta_high, n_points);
for i = 1:n_points
    u_d1(i, :) = u(theta_range(i), K1, k1);
    u_d2(i, :) = u(theta_range(i), K2, k2);
end
% figure();
% plot(theta_range, u_d1) 
% legend('d1', 'd2', 'd3')
% title('Utility Curves for 3-action numerical verification');

% figure();
% plot(theta_range, u_d2) 
% legend('d1', 'd2')
% title('Utility Curves for 2-action numerical verification');

evsi_analytical_n1 = @(n) evsi_normal_normal_analytical_n(n, N(1), mu0, n0, sigma, K1, k1, B1);
evsi_analytical_n2 = @(n) evsi_normal_normal_analytical_n(n, N(2), mu0, n0, sigma, K2, k2, B2);

cm = lines(5);
% plot evsi and cost curves
f1_values = zeros(n_points_n, 1);
f2_values = zeros(n_points_n, 1);
c_values = zeros(n_points_n, 1);
for i = 1:n_points_n
    f1_values(i) = evsi_analytical_n1(n_grid_course(i));
    f2_values(i) = evsi_analytical_n2(n_grid_course(i));
    c_values(i) = c * n_grid_course(i);
end

% obtain unconstrained optimizers for both
% only using s as an input for compatability with peicewise solver
evsi_analytical_ns1 = @(n, s) evsi_normal_normal_analytical_n(n, N(1), mu0, n0, sigma, K1, k1, B1);
evsi_analytical_ns2 = @(n, s) evsi_normal_normal_analytical_n(n, N(2), mu0, n0, sigma, K2, k2, B2);
evsi_analytical_ns12 = @(n, s) evsi_analytical_ns1(n, s) + evsi_analytical_ns2(n, s);

n_opt_uncon1 = piecewise_solver(1, evsi_analytical_ns1, c, 1e7, 20000, 1000);
n_opt_uncon2 = piecewise_solver(1, evsi_analytical_ns2, c, 1e7, 20000, 1000);
n_opt_uncon12 = piecewise_solver(1, evsi_analytical_ns12, c, 1e7, 20000, 1000);

% solve with a budget B1
n_opt_con1 = piecewise_solver(1, evsi_analytical_ns1, c, 2000000, 20000, 1000);
n_opt_con2 = piecewise_solver(1, evsi_analytical_ns2, c, 2000000, 20000, 1000);
n_opt_con12 = piecewise_solver(1, evsi_analytical_ns12, c, 2000000, 20000, 1000);

% solve with a budget B2
n_opt_con1_2 = piecewise_solver(1, evsi_analytical_ns1, c, 3000000, 20000, 1000);
n_opt_con2_2 = piecewise_solver(1, evsi_analytical_ns2, c, 3000000, 20000, 1000);
n_opt_con12_2 = piecewise_solver(1, evsi_analytical_ns12, c, 3000000, 20000, 1000);

% evsi + cost
figure
hold on
plot(n_grid_course, f1_values, 'LineWidth', 4, 'Color', cm(1, :))
plot(n_grid_course, f2_values, 'LineWidth', 4, 'Color', cm(2, :))
plot(n_grid_course, c_values, '--', 'LineWidth', 4, 'Color', cm(3, :))
text(n_u_bnd, f1_values(end)+10, '$v_1(n)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_u_bnd, f2_values(end)+10, '$v_2(n)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_u_bnd, c_values(end)+10, '$c(n)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
xlabel('Study sample size (n)', 'Interpreter', 'latex')
ylabel('Expected value (\$)', 'Interpreter', 'latex')
set(gca,'FontSize',16)
% legend('Decision Problem 1', 'Decision Problem 2', 'Cost')

% enbs and total enbs
close all
figure
hold on
plot(n_grid_course, f1_values - c_values, 'LineWidth', 4, 'Color', cm(1, :))
plot(n_grid_course, f2_values - c_values, 'LineWidth', 4, 'Color', cm(2, :))
plot(n_grid_course, f1_values + f2_values - c_values, 'LineWidth', 4, 'Color', cm(4, :))
text(n_u_bnd, f1_values(end) - c_values(end) +0.75e6, '$\mathrm{ENBS}_1(n)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_u_bnd, f2_values(end) - c_values(end) +1.2e6, '$\mathrm{ENBS}_2(n)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_u_bnd, f1_values(end) + f2_values(end) - c_values(end) +0.5e6, '$\mathrm{ENBS}(n)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')

plot(n_opt_uncon1, evsi_analytical_n1(n_opt_uncon1) - c*n_opt_uncon1, '*', 'Color', cm(5, :), 'MarkerSize', 15, 'HandleVisibility','off')
plot(n_opt_uncon2, evsi_analytical_n2(n_opt_uncon2) - c*n_opt_uncon2, '*',  'Color', cm(5, :), 'MarkerSize', 15, 'HandleVisibility','off')
plot(n_opt_uncon12, evsi_analytical_ns12(n_opt_uncon12, 1) - c*n_opt_uncon12, '*',  'Color', cm(5, :), 'MarkerSize', 15, 'HandleVisibility','off')

plot(n_opt_con1, evsi_analytical_n1(n_opt_con1) - c*n_opt_con1, '.',  'Color', cm(5, :), 'MarkerSize', 20, 'HandleVisibility','off')
plot(n_opt_con2, evsi_analytical_n2(n_opt_con2) - c*n_opt_con2, '.',  'Color', cm(5, :), 'MarkerSize', 20, 'HandleVisibility','off')
plot(n_opt_con12, evsi_analytical_ns12(n_opt_con12, 1) - c*n_opt_con12, '.',  'Color', cm(5, :), 'MarkerSize', 20, 'HandleVisibility','off')
text(n_opt_con1, evsi_analytical_n1(n_opt_con1) - c*n_opt_con1 + 0.21e6, 'B1', 'Color', cm(1, :), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(n_opt_con2, evsi_analytical_n2(n_opt_con2) - c*n_opt_con2+ 0.1e6, 'B1', 'Color', cm(2, :), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(n_opt_con12, evsi_analytical_ns12(n_opt_con12, 1) - c*n_opt_con12+ 0.1e6, 'B1', 'Color', cm(4, :), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')


plot(n_opt_con1_2, evsi_analytical_n1(n_opt_con1_2) - c*n_opt_con1_2, '.',  'Color', cm(5, :), 'MarkerSize', 20, 'HandleVisibility','off')
plot(n_opt_con2_2, evsi_analytical_n2(n_opt_con2_2) - c*n_opt_con2_2, '.',  'Color', cm(5, :), 'MarkerSize', 20, 'HandleVisibility','off')
plot(n_opt_con12_2, evsi_analytical_ns12(n_opt_con12_2, 1) - c*n_opt_con12_2, '.',  'Color', cm(5, :), 'MarkerSize', 20, 'HandleVisibility','off')
text(n_opt_con1_2, evsi_analytical_n1(n_opt_con1_2) - c*n_opt_con1_2 + 0.7e6, 'B2', 'Color', cm(1, :), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')
text(n_opt_con2_2, evsi_analytical_n2(n_opt_con2_2) - c*n_opt_con2_2 + 0.1e6, 'B2', 'Color', cm(2, :), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','center')
text(n_opt_con12_2, evsi_analytical_ns12(n_opt_con12_2, 1) - c*n_opt_con12_2 + 0.1e6, 'B2', 'Color', cm(4, :), 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right')


xlabel('Study sample size (n)', 'Interpreter', 'latex')
ylabel('Expected value (\$)', 'Interpreter', 'latex')
set(gca,'FontSize',16)






