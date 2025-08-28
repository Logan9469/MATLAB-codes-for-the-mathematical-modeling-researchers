close all; clc;
% parameter ranges
paramRanges = struct('beta', [0.016, 0.024], 'sigma', [0.0552, 0.0828], 'gamma', [0.00020, 0.00030], ...
                     'theta', [0.72, 1.08], 'phi', [0.064, 0.096], 'mu', [0.008, 0.012], ...
                     'alpha', [0.0002832, 0.0004248], 'lambda', [20.00, 30.00], 'r', [2.40, 3.60], ...
                     'u1', [0.000668, 0.001002], 'u2', [0.00796, 0.01194], 'rho', [0.00224, 0.00336]);

% % Scale sampled parameters within ±20% of their baseline values
% for i = 1:nParams
%     baseline = params.(paramNames{i});
%     scaledSamples(:, i) = baseline * (0.8 + 0.4 * samples(:, i));
% end



% Number of samples
nSamples = 1000;
paramNames = fieldnames(paramRanges);
nParams = numel(paramNames);

% Generate LHS samples
lhsSamples = lhsdesign(nSamples, nParams);

% Scale samples to parameter ranges
scaledSamples = zeros(nSamples, nParams);
for i = 1:nParams
    low = paramRanges.(paramNames{i})(1);
    high = paramRanges.(paramNames{i})(2);
    scaledSamples(:, i) = low + lhsSamples(:, i) * (high - low);
end

% Run simulations to compute R0 for each parameter set
R0Values = zeros(nSamples, 1);
for i = 1:nSamples
    R0Values(i) = calculateR0(scaledSamples(i, :));
end

% Calculate PRCC
prccValues = calculatePRCC(R0Values, scaledSamples);

% Display results
for i = 1:nParams
    fprintf('PRCC for %s: %.4f\n', paramNames{i}, prccValues(i));
end

% Scatter plot of each parameter vs R0 values
for i = 1:nParams
    figure;
    scatter(scaledSamples(:, i), R0Values, 'filled');
    xlabel(paramNames{i}, 'Interpreter', 'latex');
    ylabel('R0')
    %title(['2D Plot of ', paramNames{i}, ' vs. R0']);
end

% Perform PRCC analysis for R0 values
[prccValues, pValues] = calculatePRCC(R0Values, scaledSamples);

% Display PRCC results with p-values
for i = 1:nParams
    fprintf('PRCC for %s: %.4f (p-value: %.4f)\n', paramNames{i}, prccValues(i), pValues(i));
end

% PRCC bar chart with significance indication
figure;
bar(prccValues, 'FaceColor', [0.2, 0.6, 0.8]); % Default color for all bars
hold on;

% Highlight bars with significant p-values (e.g., p < 0.05)
significantIdx = find(pValues < 0.05);
bar(significantIdx, prccValues(significantIdx), 'FaceColor', [0.8, 0.2, 0.2]); % Red color for significant bars

% labels and title
xticks(1:nParams);
xticklabels({'$\beta$', '$\sigma$', '$\gamma$', '$\theta$', '$\phi$', '$\mu$', '$\alpha$', '$\lambda$', '$r$', '$u_1$', '$u_2$', '$\rho$'}); % Use LaTeX
set(gca, 'TickLabelInterpreter', 'latex'); % Enable LaTeX interpretation
ylabel('PRCC');
title('Partial Rank Correlation Coefficients with Significance');
legend('Not Significant', 'Significant (p < 0.05)');
hold off;



% --- Functions ---

function [prcc, pValues] = calculatePRCC(outputs, samples)
    nParams = size(samples, 2);
    prcc = zeros(nParams, 1);
    pValues = zeros(nParams, 1);
    for i = 1:nParams
        % Rank data
        ranksSamples = tiedrank(samples(:, i));
        ranksOutputs = tiedrank(outputs);
        % Compute correlation and p-value
        [prcc(i), pValues(i)] = corr(ranksSamples, ranksOutputs, 'Type', 'Spearman');
    end
end

% Simulate model over time for each parameter set
timeSteps = 1000; % Define time steps
Article2R0TimeSeries = zeros(nSamples, timeSteps);
for i = 1:nSamples
    Article2R0TimeSeries(i, :) = calculateR0OverTime(scaledSamples(i, :), timeSteps);
end

% number of time steps
timeSteps = 1000;

% Simulate model over time for each parameter set
Article2R0TimeSeries = zeros(nSamples, timeSteps);
for i = 1:nSamples
    Article2R0TimeSeries(i, :) = calculateR0OverTime(scaledSamples(i, :), timeSteps);
end

% Calculate time-varying PRCC
prccTimeSeries = zeros(nParams, timeSteps);
pValuesTimeSeries = zeros(nParams, timeSteps);
for t = 1:timeSteps
    [prccTimeSeries(:, t), pValuesTimeSeries(:, t)] = calculatePRCC(Article2R0TimeSeries(:, t), scaledSamples);
end

% Plot PRCC over time for each parameter
figure;
for i = 1:nParams
    plot(1:timeSteps, prccTimeSeries(i, :), 'DisplayName', paramNames{i}); % Keep auto labels for now
    hold on;
end
xlabel('Time');
ylabel('PRCC');
title('Time-varying PRCC of Each Parameter with R_0', 'Interpreter', 'latex');
legend({'$\beta$', '$\sigma$', '$\gamma$', '$\theta$', '$\phi$', '$\mu$', '$\alpha$', '$\lambda$', '$r$', '$u_1$', '$u_2$', '$\rho$'}, ...
    'Location', 'best', 'Interpreter', 'latex');% Use LaTeX for symbols
grid minor;
hold off;


function R0 = calculateR0(params)
    beta = params(1);
    sigma = params(2);
    gamma = params(3);
    theta = params(4);
    phi   = params(5);
    mu    = params(6);
    alpha = params(7);
    lambda = params(8);
    r     = params(9);
    u1     = params(10);
    u2     = params(11);
    rho     = params(12);
    % Basic reproduction number for SACR model
     R0 = (lambda * gamma) / ((mu + u1) * (sigma + theta + mu)) + ...
     (2 * beta * theta) / (((rho + (r * u2) / (1 + phi^2)) + mu - lambda * alpha) * (sigma + theta + mu));
end
