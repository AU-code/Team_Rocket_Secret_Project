%% Final project - part 3

% Seasonal Influenza
beta1 = 0.3; % Transmission Rate
gamma1 = 0.1; % Recovery Rate

% COVID
beta2 = 1; % Transmission Rate
gamma2 = 0.1; % Recovery Rate

% Measles
beta3 = 2; % Transmission Rate
gamma3 = 0.2; % Recovery Rate

% Initial values
h = 1; % time step, days
N = 1000; % Total Population, constant
t = 0:h:100; % total simulation time, 100 days
S0 = 990; % Initial Condition for S(t)
I0 = 10; % Initial Condition for I(t)
R0 = 0; % Initial Condition for R(t)
y0 = [S0 I0 R0]; % initial conditions vector
T = 100; % total simulation time
days_30 = 30; % days for linear least squares (30)
days_10 = 10; % days for linear least squares (10)


% Scenarios
scenarios = [0.3, 0.1; 1.0, 0.1; 2.0, 0.2]; % [β, γ] for Influenza, COVID, Measles
scenario_names = {'Seasonal Influenza', 'COVID', 'Measles'};

% Initialize storage for results
results = cell(size(scenarios, 1), 4); % Columns: Scenario, I0_30, β_30, I0_10, β_10

% Run simulations for each scenario
for s = 1:size(scenarios, 1)
    beta = scenarios(s, 1);
    gamma = scenarios(s, 2);

    % Initialize variables
    S = zeros(T+1, 1);
    I = zeros(T+1, 1);
    R = zeros(T+1, 1);
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;

    % Simulate the SIR model
    for t = 1:T
        dS = -beta * S(t) * I(t) / N;
        dI = beta * S(t) * I(t) / N - gamma * I(t);
        dR = gamma * I(t);
        
        S(t+1) = S(t) + h * dS;
        I(t+1) = I(t) + h * dI;
        R(t+1) = R(t) + h * dR;
    end

    % Extract "true" I(t) data for 30 days
    t30 = (1:days_30)';
    I_data_30 = I(1:days_30);

    % Apply linear least squares using Eq. (8)
    ln_I_data_30 = log(I_data_30); % Take natural log of I(t)
    X_30 = [ones(length(t30), 1), t30]; % Design matrix for [ln I(0), k]
    params_30 = X_30 \ ln_I_data_30; % Linear regression
    ln_I0_30 = params_30(1);
    k_30 = params_30(2);

    % Estimate I(0) and beta for 30 days
    I0_est_30 = exp(ln_I0_30);
    beta_est_30 = (k_30 + gamma) * N / S0;

    % Repeat for 10 days of data
    t10 = (1:days_10)';
    I_data_10 = I(1:days_10);
    ln_I_data_10 = log(I_data_10); % Take natural log of I(t)
    X_10 = [ones(length(t10), 1), t10]; % Design matrix
    params_10 = X_10 \ ln_I_data_10; % Linear regression
    ln_I0_10 = params_10(1);
    k_10 = params_10(2);

    % Estimate I(0) and beta for 10 days
    I0_est_10 = exp(ln_I0_10);
    beta_est_10 = (k_10 + gamma) * N / S0;

    % Save results
    results{s, 1} = scenario_names{s};
    results{s, 2} = I0_est_30;
    results{s, 3} = beta_est_30;
    results{s, 4} = I0_est_10;
    results{s, 5} = beta_est_10;
end

% Display results
fprintf('Results:\n');
fprintf('%-20s %-10s %-10s %-10s %-10s\n', 'Scenario', 'I0_30', 'β_30', 'I0_10', 'β_10');
for s = 1:size(scenarios, 1)
    fprintf('%-20s %-10.2f %-10.2f %-10.2f %-10.2f\n', results{s, 1}, results{s, 2}, results{s, 3}, results{s, 4}, results{s, 5});
end
