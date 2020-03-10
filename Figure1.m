% Author: Michael Fairley
% Date: December 11, 2019
% MATLAB 2019b Update 2

% use a logit for envsi and subtract cost
K = 1000;
A = 10;
k = 0.08;
b = 60;
cf = 0;

logit = @(t) K ./ (1 + A .* exp(-k .* (t - b)));
% inflection point
ti = log(A) / k + 60;
evsi_ti = logit(ti);

% set cost to make enbs 0 exactly at inflection point
% evsi_ti - c * ti = 0
c = (evsi_ti - cf) / ti;
c_func = @(t) c * t + (t > 0) * cf;
enbs = @(t) logit(t) - c_func(t);
enbs_neg = @(t) - enbs(t);
% maximizing point
tstar = fmincon(enbs_neg, 125, [], [], [], [], 0);

figure
hold on
n_ub = 250;
axis_low = -600;
axis_high = 1600;
line_width = 3;
fplot(logit, [0, n_ub], 'k', 'LineWidth', line_width)
fplot(enbs, [0, n_ub], 'k', 'LineWidth', line_width)
fplot(c_func, [0, n_ub], '--k', 'LineWidth', line_width)
plot(tstar, enbs(tstar), '*k', 'MarkerSize', 15)
yline(0)
% dotted lines
% vertical
plot([ti, ti], [0, evsi_ti], ':k', 'LineWidth', line_width); 
plot([0, ti], [evsi_ti, evsi_ti], ':k', 'LineWidth', line_width); 
axis([0, n_ub, axis_low, axis_high])
% labels - lines
text(n_ub, c_func(n_ub)+2, '$c_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_ub, logit(n_ub)+10, '$N_s v_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
text(n_ub, enbs(n_ub)+230, '$\mathrm{ENBS}_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
% labels - points
text(ti, 0, '$\underline{n_s}$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','top','HorizontalAlignment','left')
text(0, evsi_ti, '$c_s^f$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','middle','HorizontalAlignment','right')
xlabel('Study sample size $(n_s)$', 'Interpreter', 'latex')
ylabel('Expected value (\$)', 'Interpreter', 'latex')
set(gca,'FontSize',16)

% second panel with separated cost functions
% new continuous cost function
% cc_func = @(t) c_func(t) - c_func(ti);
% figure
% hold on
% fplot(logit, [0, n_ub], 'k', 'LineWidth', line_width)
% fplot(enbs, [0, n_ub], 'k', 'LineWidth', line_width)
% yline(c_func(ti), '--k', 'LineWidth', line_width)
% fplot(cc_func, [0, n_ub], '--k', 'LineWidth', line_width)
% plot(tstar, enbs(tstar), '*k', 'MarkerSize', 15)
% yline(0)
% % dotted lines
% % vertical
% axis([0, n_ub, axis_low, axis_high])
% % labels - lines
% % text(n_ub, c_func(n_ub)+2, '$c_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
% text(n_ub, logit(n_ub)+10, '$N_s v_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
% text(n_ub, enbs(n_ub)+230, '$\mathrm{ENBS}_s(n_s)$', 'FontSize',16, 'VerticalAlignment','bottom','HorizontalAlignment','right', 'Interpreter', 'latex')
% text(n_ub, c_func(ti)+60, '$c_s^f$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','middle','HorizontalAlignment','right')
% text(n_ub, cc_func(n_ub)+20, '$c_s^c(n_s)$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','middle','HorizontalAlignment','right')
% % labels - points
% text(ti, 0, '$\underline{n_s}$', 'Interpreter', 'latex', 'FontSize',16, 'VerticalAlignment','top','HorizontalAlignment','left')
% xlabel('Study sample size $(n_s)$', 'Interpreter', 'latex')
% ylabel('Expected value (\$)', 'Interpreter', 'latex')
% set(gca,'FontSize',16)



