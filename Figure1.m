% Author: Michael Fairley
% Date: December 11, 2019
% MATLAB 2019b Update 2

% use a logit for envsi and subtract cost
K = 1000;
A = 10;
k = 0.08;

logit = @(t) K ./ (1 + A .* exp(-k .* (t - 60)));
% inflection point
ti = log(A) / k + 60;
evsi_ti = logit(ti);

% set cost to make enbs 0 exactly at inflection point
% evsi_ti - c * ti = 0
c = evsi_ti / ti;
c_func = @(t) c * t;
enbs = @(t) logit(t) - c * t;
figure
hold on
n_ub = 250;
fplot(logit, [0, n_ub], 'k', 'LineWidth', 4)
fplot(enbs, [0, n_ub], 'k', 'LineWidth', 4)
fplot(c_func, [0, n_ub], '--k', 'LineWidth', 4)
yline(0)
% dotted lines
% vertical
plot([ti, ti], [0, evsi_ti], ':k', 'LineWidth', 4); 
plot([0, ti], [evsi_ti, evsi_ti], ':k', 'LineWidth', 4); 
axis([0, n_ub, -400, 1600])
% labels - lines
text(n_ub, c_func(n_ub)+10, '$c_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_ub, logit(n_ub)+10, '$v_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_ub, enbs(n_ub)+230, '$\mathrm{ENBS}_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
% labels - points
text(ti, 0, '$\underline{n_s}$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','top','HorizontalAlignment','left')
text(0, evsi_ti, '$c_s^f$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','middle','HorizontalAlignment','right')
xlabel('Study sample size $(n_s)$', 'Interpreter', 'latex')
ylabel('Expected value (\$)', 'Interpreter', 'latex')
set(gca,'FontSize',16)