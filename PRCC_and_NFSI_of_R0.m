%% Sensitivity Analysis of R0 (Normalized and PRCC)
clc; clear; close all;

%% Step 1: Define symbols
syms A beta1 beta2 phi nu alpha theta delta mu psi a tau rho epsilon

%% Step 2: Define R0 expression
numerator = nu*A*(mu + delta + phi)*(tau*beta1 + beta2*(epsilon + theta - rho));
denominator = (mu + alpha + a)*(epsilon + theta - rho)*((mu + delta + phi)*(mu + psi) - phi*psi);
R0 = numerator / denominator;

%% Step 3: Parameter values
param_syms = [ A, beta1, beta2, phi, nu, alpha, theta, delta, mu, psi, a, tau, rho, epsilon ];
param_vals = [ 0.9, 0.0414, 0.00115, 0.0051, 0.01004, 0.0018, 0.9, 0.89, ...
               0.004563, 0.5, 0.010714, 0.9, 0.33, 0.1724 ];
param_names = {'$A$','$\beta_1$','$\beta_2$','$\phi$','$\nu$','$\alpha$','$\theta$', ...
               '$\delta$','$\mu$','$\psi$','$a$','$\tau$','$\rho$','$\epsilon$'};

%% Step 4: Compute R0 value numerically
R0_val = double(subs(R0, param_syms, param_vals));

%% Step 5: Compute Normalized Sensitivity Indices
SI_vals = zeros(1, length(param_syms));
for i = 1:length(param_syms)
    dR0_dp = diff(R0, param_syms(i));
    SI_expr = simplify(dR0_dp * param_syms(i) / R0);
    SI_vals(i) = double(subs(SI_expr, param_syms, param_vals));
end

%% Step 6: Bar plot - Vertical for Normalized Sensitivity Index
figure;
bar_colors = repmat([1 0 0], length(SI_vals), 1);
bar_colors(SI_vals >= 0, :) = repmat([0 0.6 0], sum(SI_vals >= 0), 1);
b = bar(SI_vals, 'FaceColor', 'flat', 'BarWidth', 0.6);
b.CData = bar_colors;
for i = 1:length(SI_vals)
    value = SI_vals(i);
    ypos = value + 0.02 * sign(value);
    if value >= 0
        align = 'bottom';
    else
        align = 'top';
    end
    text(i, ypos, sprintf('%.2f', value), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', align, 'FontSize', 10);
end
set(gca, 'XTick', 1:numel(param_names), 'XTickLabel', param_names, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 12);
ylabel('Sensitivity Index', 'Interpreter', 'latex');
title('Normalized Sensitivity Indices of $\mathcal{R}_0$', 'Interpreter', 'latex');
grid minor; yline(0, '--k');

%% Step 7: PRCC Analysis with P-values
N = 1000;  % samples
lb = 0.8 * param_vals;
ub = 1.2 * param_vals;
lhs = lhsdesign(N, length(param_vals));
samples = lhs .* (ub - lb) + lb;

% Evaluate R0
R0_vals = zeros(N,1);
for i = 1:N
    p = num2cell(samples(i,:));
    [A, beta1, beta2, phi, nu, alpha, theta, delta, mu, psi, a, tau, rho, epsilon] = deal(p{:});
    num = nu*A*(mu + delta + phi)*(tau*beta1 + beta2*(epsilon + theta - rho));
    den = (mu + alpha + a)*(epsilon + theta - rho)*((mu + delta + phi)*(mu + psi) - phi*psi)*(mu + nu);
    R0_vals(i) = num / den;
end

% PRCC calculation
ranked_params = tiedrank(samples);
ranked_R0 = tiedrank(R0_vals);
PRCC_vals = zeros(1, length(param_names));
p_vals = zeros(1, length(param_names));
for i = 1:length(param_names)
    Xi = ranked_params(:, i);
    X_rest = ranked_params;
    X_rest(:, i) = [];
    r_Xi = Xi - X_rest * (X_rest \ Xi);
    r_Y  = ranked_R0 - X_rest * (X_rest \ ranked_R0);
    [R, P] = corrcoef(r_Xi, r_Y);
    PRCC_vals(i) = R(1,2);
    p_vals(i) = P(1,2);
end

%% Step 8: PRCC Bar Plot with Color
figure;
bar_colors_prcc = repmat([1 0 0], length(PRCC_vals), 1);
bar_colors_prcc(PRCC_vals >= 0, :) = repmat([0 0.6 0], sum(PRCC_vals >= 0), 1);
b2 = bar(PRCC_vals, 'FaceColor', 'flat', 'BarWidth', 0.6);
b2.CData = bar_colors_prcc;
for i = 1:length(PRCC_vals)
    ypos = PRCC_vals(i) + 0.02 * sign(PRCC_vals(i));
    if PRCC_vals(i) >= 0
        align = 'bottom';
    else
        align = 'top';
    end
    text(i, ypos, sprintf('%.2f', PRCC_vals(i)), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', align, 'FontSize', 10);
end
set(gca, 'XTick', 1:numel(param_names), 'XTickLabel', param_names, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 12);
ylabel('PRCC', 'Interpreter', 'latex');
title('Partial Rank Correlation Coefficients (PRCC) of $\mathcal{R}_0$', 'Interpreter', 'latex');
grid minor; yline(0, '--k');
%print('PRCC Bar Graph', '-dpng', '-r800');

%% Step 10: Scatter Plots of Each Parameter vs R0
for i = 1:length(param_syms)
    figure;
    scatter(samples(:, i), R0_vals, 20, 'filled');
    xlabel(param_names{i}, 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\mathcal{R}_0$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['2D Scatter Plot of ', param_names{i}, ' vs. $\mathcal{R}_0$'], ...
          'Interpreter', 'latex', 'FontSize', 13);
    grid minor;
end


%% Step 9: Export LaTeX Table for PRCC
fprintf('\n%% LaTeX Table Format (PRCC with p-values)\n');
for i = 1:length(param_names)
    cleanname = erase(param_names{i}, {'$', '\'}); % remove LaTeX for raw table
    fprintf('%s & %.4f & %.4g \\\n', cleanname, PRCC_vals(i), p_vals(i));
end
