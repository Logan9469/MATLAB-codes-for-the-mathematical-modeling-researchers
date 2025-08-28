function R0TimeSeries = calculateR0OverTime(params, timeSteps)
    % Parameters
    beta = params(1);
    sigma = params(2);
    gamma = params(3);
    theta = params(4);
    tau   = params(5);
    mu    = params(6);
    alpha = params(7);
    lambda = params(8);
    r     = params(9);
    u     = params(10);

    % Initialize R0 time series
    R0TimeSeries = zeros(1, timeSteps);

    % Loop over time steps to simulate R0 over time
    for t = 1:timeSteps
        % Example: Suppose tau changes over time as an example of time-dependence
        tau_t = tau * (1 + 0.1 * sin(2 * pi * t / timeSteps)); % Example oscillation
        % Calculate R0 with time-varying tau
        R0TimeSeries(t) = 2 * beta * (theta + tau_t + mu - alpha * lambda + sigma * gamma) ...
                          / ((sigma + mu + r * u) * (theta + tau_t + mu - alpha * lambda));
    end
end
