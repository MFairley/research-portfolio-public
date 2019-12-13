% Author: Michael Fairley
% Date: December 11, 2019
% MATLAB 2019b Update 2
n01 = 5;
n02 = 5;
rho0 = 0.0;
rho05 = 0.9;
rho09 = 1.0;

nu_corr = @(n1, n2, rho) (n1 * n02 + n1*n2*(1-rho^2) + n2*n01*rho^2) / ...
    (n1*n02 + n2*n01 + n1*n2*(1-rho^2) + n01*n02);

% generate data - rho 0
points = 1000;
x = linspace(0, 20, points); % n1
y = linspace(0, 25, points); % n2
[X, Y] = meshgrid(x, y);

Z = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        Z(i, j) = sqrt(nu_corr(X(i, j), Y(i, j), rho0));
    end
end

contourf(X, Y, Z, 'ShowText','on') % mi
xlabel('n_1');
ylabel('n_2');
caxis manual
title(sprintf('$$\\rho = %.2f$$', rho0), 'interpreter', 'latex')
caxis([0 1.00]);
axis square
set(gca,'FontSize',16)

% generate data - rho 0.5
points = 1000;
x = linspace(0, 20, points); % n1
y = linspace(0, 20, points); % n2
[X, Y] = meshgrid(x, y);

Z = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        Z(i, j) = sqrt(nu_corr(X(i, j), Y(i, j), rho05));
    end
end

figure
contourf(X, Y, Z, 'ShowText','on') % mi
xlabel('n_1');
ylabel('n_2');
caxis manual
title(sprintf('$$\\rho = %.2f$$', rho05), 'interpreter', 'latex')
caxis([0 1.00]);
axis square
set(gca,'FontSize',16)


% generate data - rho 0.5
points = 1000;
x = linspace(0, 20, points); % n1
y = linspace(0, 20, points); % n2
[X, Y] = meshgrid(x, y);

Z = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        Z(i, j) = sqrt(nu_corr(X(i, j), Y(i, j), rho09));
    end
end

figure
contourf(X, Y, Z, 'ShowText','on') % mi
xlabel('n_1');
ylabel('n_2');
title(sprintf('$$\\rho = %.2f$$', rho09), 'interpreter', 'latex')
caxis manual
caxis([0 1.00]);
axis square
set(gca,'FontSize',16)
