%% Horizontal Sensitivity Analysis of R0
clc; clear; close all;

% Step 1: Define symbols
syms A beta1 beta2 phi nu alpha theta delta mu psi a tau rho epsilon

% Step 2: Define R0 expression
numerator = nu*A*(mu + delta + phi)*(tau*beta1 + beta2*(epsilon + theta - rho));
denominator = (mu + alpha + a)*(epsilon + theta - rho)*((mu + delta + phi)*(mu + psi) - phi*psi);
R0 = numerator / denominator;

% Step 3: Parameter values
param_syms = [ A, beta1, beta2, phi, nu, alpha, theta, delta, mu, psi, a, tau, rho, epsilon ];
param_vals = [ 0.9, 0.0414, 0.00115, 0.0051, 0.01004, 0.0018, 0.9, 0.89, ...
               0.004563, 0.5, 0.010714, 0.9, 0.33, 0.1724 ];

% LaTeX-style labels for plot
param_latex = {'$A$', '$\beta_1$', '$\beta_2$', '$\phi$', '$\nu$', '$\alpha$', ...
               '$\theta$', '$\delta$', '$\mu$', '$\psi$', '$a$', '$\tau$', ...
               '$\rho$', '$\epsilon$'};

% Step 4: Compute R0 value numerically
R0_val = double(subs(R0, param_syms, param_vals));

% Step 5: Compute Sensitivity Indices
SI_vals = zeros(1, length(param_syms));
for i = 1:length(param_syms)
    dR0_dp = diff(R0, param_syms(i));  % partial derivative
    SI_expr = simplify(dR0_dp * param_syms(i) / R0);  % normalized index
    SI_vals(i) = double(subs(SI_expr, param_syms, param_vals));
    fprintf('Sensitivity index for %-7s = %+0.4f\n', char(param_syms(i)), SI_vals(i));
end

% Step 6: Assign bar colors
bar_colors = repmat([1 0 0], length(SI_vals), 1);  % red (default: negative)
bar_colors(SI_vals >= 0, :) = repmat([0 0.6 0], sum(SI_vals >= 0), 1);  % green if positive

% Step 7: Plot Horizontal Bar Chart
figure;
b = barh(SI_vals, 'FaceColor', 'flat');
b.CData = bar_colors;

% Add numeric labels to bars
for i = 1:length(SI_vals)
    value = SI_vals(i);
    xpos = value + 0.01 * sign(value);  % offset from bar
    if value >= 0
        align = 'left';
    else
        align = 'right';
    end
    text(xpos, i, sprintf('%.2f', value), 'VerticalAlignment', 'middle', ...
         'HorizontalAlignment', align, 'FontSize', 10);
end

% Formatting and axis settings
set(gca, 'YTick', 1:numel(param_latex), 'YTickLabel', param_latex, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 12);
xlabel('Sensitivity Index', 'Interpreter', 'latex', 'FontSize', 14);
%title('Normalized Sensitivity Indices of $R_0$', 'Interpreter', 'latex', 'FontSize', 15);
grid minor;
xline(0, '--k', 'LineWidth', 1);
%print('Sensitivity Indices', '-dpng', '-r800');


%% Vertical Bar graph of R0

%% Vertical Sensitivity Analysis of R0
clc; clear; close all;

% Step 1: Define symbols
syms A beta1 beta2 phi nu alpha theta delta mu psi a tau rho epsilon

% Step 2: Define R0 expression
numerator = nu*A*(mu + delta + phi)*(tau*beta1 + beta2*(epsilon + theta - rho));
denominator = (mu + alpha + a)*(epsilon + theta - rho)*((mu + delta + phi)*(mu + psi) - phi*psi);
R0 = numerator / denominator;

% Step 3: Parameter values
param_syms = [ A, beta1, beta2, phi, nu, alpha, theta, delta, mu, psi, a, tau, rho, epsilon ];
param_vals = [ 0.9, 0.0414, 0.00115, 0.0051, 0.01004, 0.0018, 0.9, 0.89, ...
               0.004563, 0.5, 0.010714, 0.9, 0.33, 0.1724 ];

% LaTeX-style labels for plot
param_latex = {'$A$', '$\beta_1$', '$\beta_2$', '$\phi$', '$\nu$', '$\alpha$', ...
               '$\theta$', '$\delta$', '$\mu$', '$\psi$', '$a$', '$\tau$', ...
               '$\rho$', '$\epsilon$'};

% Step 4: Compute R0 value numerically
R0_val = double(subs(R0, param_syms, param_vals));

% Step 5: Compute Sensitivity Indices
SI_vals = zeros(1, length(param_syms));
for i = 1:length(param_syms)
    dR0_dp = diff(R0, param_syms(i));  % partial derivative
    SI_expr = simplify(dR0_dp * param_syms(i) / R0);  % normalized index
    SI_vals(i) = double(subs(SI_expr, param_syms, param_vals));
    fprintf('Sensitivity index for %-7s = %+0.4f\n', char(param_syms(i)), SI_vals(i));
end

% Step 6: Assign bar colors
bar_colors = repmat([1 0 0], length(SI_vals), 1);  % red (default: negative)
bar_colors(SI_vals >= 0, :) = repmat([0 0.6 0], sum(SI_vals >= 0), 1);  % green if positive

% Step 7: Plot Vertical Bar Chart
figure;
b = bar(SI_vals, 'FaceColor', 'flat', 'BarWidth', 0.6);  % narrower bars
b.CData = bar_colors;

% Add numeric labels to bars
for i = 1:length(SI_vals)
    value = SI_vals(i);
    ypos = value + 0.02 * sign(value);  % larger offset from bar
    if value >= 0
        align = 'bottom';
    else
        align = 'top';
    end
    text(i, ypos, sprintf('%.2f', value), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', align, 'FontSize', 10);
end

% Formatting and axis settings
set(gca, 'XTick', 1:numel(param_latex), 'XTickLabel', param_latex, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 12);
ylabel('Sensitivity Index', 'Interpreter', 'latex', 'FontSize', 14);
%title('Normalized Sensitivity Indices of $R_0$', 'Interpreter', 'latex', 'FontSize', 15);
grid minor;
yline(0, '--k', 'LineWidth', 1);
%print('Sensitivity Indices', '-dpng', '-r800');

