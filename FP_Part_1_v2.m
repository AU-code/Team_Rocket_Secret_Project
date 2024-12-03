%% Final Project - Part 1

% Seasonal Influenza
beta_influenza = 0.3; % Transmission Rate
gamma_influenza = 0.1; % Recovery Rate

% COVID
beta_covid = 1; % Transmission Rate
gamma_covid = 0.1; % Recovery Rate

% Measles
beta_measles = 2; % Transmission Rate
gamma_measles = 0.2; % Recovery Rate

h = 1; % Time Step (1 day)
N = 1000; % Total Population (constant)
t = 0:h:100; % Time Span (100 days)
S0 = 990; % Initial Susceptible Population
I0 = 10;  % Initial Infected Population
R0 = 0;   % Initial Recovered Population


% Solve SIR Model for each disease
[S_influenza, I_influenza, R_influenza] = RK4_SIR(beta_influenza, gamma_influenza, t, N);

[S_covid, I_covid, R_covid] = RK4_SIR(beta_covid, gamma_covid, t, N);

[S_measles, I_measles, R_measles] = RK4_SIR(beta_measles, gamma_measles, t, N);

%% Plotting it ALL

% Seasonal Influenza
figure(1);
plot(t,S_influenza,'r-');
hold on
plot(t,I_influenza,'g-')
plot(t,R_influenza,'b')
hold off
xlabel('time (days)');
ylabel('Population');
title('Seasonal Influenza Simulation')
legend('S(t)','I(t)','R(t)');

% COVID
figure(2);
plot(t,S_covid,'r-');
hold on
plot(t,I_covid,'g-')
plot(t,R_covid,'b')
hold off
xlabel('time (days)');
ylabel('Population');
title('COVID Simulation')
legend('S(t)','I(t)','R(t)');

% Measles
figure(3);
plot(t,S_measles,'r-');
hold on
plot(t,I_measles,'g-')
plot(t,R_measles,'b')
hold off
xlabel('time (days)');
ylabel('Population');
title('Measles Simulation')
legend('S(t)','I(t)','R(t)');



%% Sanity Checks . . .
covid_check1 = S_covid(4) + I_covid(4) + R_covid(4)
covid_check2 = S_covid(20) + I_covid(20) + R_covid(20)
covid_check3 = S_covid(69) + I_covid(69) + R_covid(69)

flu_check1 = S_influenza(4) + I_influenza(4) + R_influenza(4)
flu_check2 = S_influenza(20) + I_influenza(20) + R_influenza(20)
flu_check3 = S_influenza(69) + I_influenza(69) + R_influenza(69)

measles_check1 = S_measles(4) + I_measles(4) + R_measles(4)
measles_check2 = S_measles(20) + I_measles(20) + R_measles(20)
measles_check3 = S_measles(69) + I_measles(69) + R_measles(69)



%% RK4 Solver for the SIR Model
function [S, I, R] = RK4_SIR(beta, gamma, t, N)

    h = 1; % time step (1 day)
    n = length(t); % n time steps
    S = zeros(n, 1); % Initialize S vector
    I = zeros(n, 1); % Initialize I vector
    R = zeros(n, 1); % Initialize R bector

    S(1) = 990; % initial susceptible
    I(1) = 10; % initial infected
    R(1) = 0; % initial recovered

    % Solve the SIR equations (ODEs) using the RK4 method
    for i = 1:(n-1)
        
        k1_S = -beta*S(i)*I(i) / N;
        k1_I = ( beta*S(i)*I(i) / N ) - ( gamma * I(i) );
        k1_R = gamma*I(i);

        k2_S = -beta * ( S(i) + h * k1_S / 2 ) * ( I(i) + h * k1_I / 2 ) / N;
        k2_I = beta * ( S(i) + (h * k1_S / 2) ) * ( I(i) + (h * k1_I / 2) ) / N - ( gamma * (I(i) + (h * k1_I / 2)) );
        k2_R = gamma * ( I(i) + h * k1_I / 2 );

        k3_S = -beta * (S(i) + (h * k2_S / 2)) * (I(i) + (h * k2_I / 2)) / N;
        k3_I = beta * ( (S(i) + (h * k2_S / 2)) * (I(i) + (h * k2_I / 2)) / N ) - ( gamma * (I(i) + (h * k2_I / 2)) );
        k3_R = gamma * (I(i) + h * k2_I / 2);

        k4_S = -beta * (S(i) + h * k3_S) * (I(i) + h * k3_I) / N;
        k4_I = beta * ( (S(i) + h * k3_S) * (I(i) + h * k3_I) / N ) - ( gamma * (I(i) + h * k3_I) );
        k4_R = gamma * (I(i) + h * k3_I);

        % updated values
        S(i+1) = S(i) + (h/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
        I(i+1) = I(i) + (h/6) * (k1_I + 2*k2_I + 2*k3_I + k4_I);
        R(i+1) = R(i) + (h/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
    end
end