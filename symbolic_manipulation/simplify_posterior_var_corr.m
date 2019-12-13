% Performs the simplification used in Appendix C of the paper
syms n1 n01 n2 n02 rho sigma1 sigma2

expr = n1*(n01+n1)*(n02+n2-n2*rho^2)^2;

expr2 = n01^2*n2*rho^2*(n02+n2);

expr3 = 2*n01*rho^2*n2*(n1*(n2+n02)-rho^2*n1*n2);

a = (n1+n01)*(n2+n02)-n1*n2*rho^2;


pretty(simplify((expr + expr2 + expr3) / a^2, 'Steps', 10000, 'IgnoreAnalyticConstraints',true))

% simplify the terms in variance
v_expr1 = (((n1*(n2+n02)-rho^2*n1*n2))^2 * (sigma1^2/n1 + sigma1^2/n01));

v_expr2 = ((rho*sigma1*sqrt(n02)*n2*sqrt(n01)/sigma2)^2 * (sigma2^2/n2 + sigma2^2/n02));

v_expr3 = (n1*(n2+n02)-rho^2*n1*n2) * (rho*sigma1*sqrt(n02)*n2*sqrt(n01))/sigma2 * (rho * sigma1 * sigma2 / sqrt(n01 * n02));

pretty(simplify(v_expr1, 'Steps', 1000))

pretty(simplify(v_expr2, 'Steps', 1000))

pretty(simplify(v_expr3, 'Steps', 1000))
