%% Final Project - Part 3

% Initial conditions
h = 1; % time step, days
N = 1000; % Total Population, constant
S0 = 990; % Initial Condition for S(t)
I0 = 10; % Initial Condition for I(t)
R0 = 0; % Initial Condition for R(t)
T = 100; % total simulation time
days_30 = 30; % days for linear least squares (30 days)
days_10 = 10; % days for linear least squares (10 days)

% Disease
disease = [0.3, 0.1; 1.0, 0.1; 2.0, 0.2]; % [beta, gamma]
disease_names = {'Seasonal Influenza', 'COVID-19', 'Measles'};

% Loop through each disease
for i = 1:length(disease_names)
    beta = disease(i, 1);
    gamma = disease(i, 2);

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

    % True I(t) data for 30 days
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

    % Display results
    fprintf('Desease: %s\n', disease_names{i});
    fprintf('True beta = %.2f, True gamma = %.2f\n', beta, gamma);
    fprintf('--- Using 30 days of data ---\n');
    fprintf('Estimated I(0): %.2f, Estimated beta: %.2f\n', I0_est_30, beta_est_30);
    fprintf('--- Using 10 days of data ---\n');
    fprintf('Estimated I(0): %.2f, Estimated beta: %.2f\n\n', I0_est_10, beta_est_10);
end
