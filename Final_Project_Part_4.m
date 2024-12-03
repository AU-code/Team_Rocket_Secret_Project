%% PART VI: Fourier Analysis
clc;
clear;

%% Part 1
B0 = 0.3; % Base transmission rate
A = 5; % Amplitude of fluctuation
w = 2 * pi; % Angular frequency
Beta_t = @(t) B0 * (1 + A * sin(w * t)); % Time-dependent beta

%% Part 2
gammaP4 = 0.1; % Recovery rate
h = 0.1; % Time step, days
N = 1000; % Total population (constant)
t = 0:h:30; % Time vector (30 days)

% Solve SIR Model with fluctuating beta
[S, I, R] = RK4_SIR(Beta_t, gammaP4, t, N);

%% Part 3
% Plot results
figure;
plot(t,S,'b',t,I,'r',t,R,'g') %Graph of SIR Model
grid on
title('SIR Model Plot for Periodic Transmission Rate (\beta(t))')
xlabel('Time (days)')
ylabel('Population')
legend('Susceptible Pop','Infected Pop','Recovered Pop')

%Yes, there are periodic fluctuations in the signals as a result of the
%periodicity of Beta.

%% Parts 4 and 5
% Spectral Analysis on Infected Population I(t)
N_spectrum = length(I);
Fs = 1 / h; % Sampling frequency (samples per day)
spectrum1 = fft(R);
spectrum2 = fft(S);
spectrum3 = fft(I); % Compute FFT of I(t)
P2 = abs(spectrum3 / N_spectrum); % Normalize and get the two-sided spectrum
P1 = P2(1:floor(N_spectrum/2) + 1); % Get the one-sided spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Scale appropriately for one-sided spectrum

% Frequency range
f = Fs * (0:(N_spectrum/2)) / N_spectrum; % Frequency vector

% Plot the spectrum
figure;
semilogy(f, P1)
title('Spectrum of Infected Population I(t)');
xlabel('Frequency (1/day)');
ylabel('Magnitude');
grid on;

% Yes it makes sense physically because the spectrum shows the periodic
% variation that happened when I changed Beta. 

%% Part 6

B0 = 0.3;%
A = 5;%Amplitude
w = 2*pi*(100/365);%frequency 
Beta_t = @(t) B0*(1+A*sin(w*t)); %Transmition rate fluctuation

gammaP4 = 0.1; % Recovery Rate

h = 0.1; % time step, days
N = 1000; % Total Population, constant
t = 0:h:30; % total simulation time, 30 days
S0 = 990; % Initial Condition for S(t)
I0 = 10; % Initial Condition for I(t)
R0 = 0; % Initial Condition for R(t)

% SOLVE SIR MODEL USING RUNGE KUTTA 4th ORDER METHOD

[S, I, R] = RK4_SIR(Beta_t, gammaP4, t, N);

figure;
plot(t,S,'b',t,I,'r',t,R,'g') %Graph of SIR Model
grid on
title('SIR Model Plot for Periodic Transmission Rate (\beta(t))')
xlabel('Time (days)')
ylabel('Population')
legend('Susceptible Pop','Infected Pop','Recovered Pop')

N_spectrum = length(I);
Fs = 1 / h; % Sampling frequency (samples per day)
spectrum1 = fft(R);
spectrum2 = fft(S);
spectrum3 = fft(I); % Compute FFT of I(t)
P2 = abs(spectrum3 / N_spectrum); % Normalize and get the two-sided spectrum
P1 = P2(1:floor(N_spectrum/2) + 1); % Get the one-sided spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Scale appropriately for one-sided spectrum

% Frequency range
f = Fs * (0:(N_spectrum/2)) / N_spectrum; % Frequency vector

% Plot the spectrum
figure;
semilogy(f, P1)
title('Spectrum of Infected Population I(t)');
xlabel('Frequency (1/day)');
ylabel('Magnitude');
grid on;

% The peak frequency shifts to a lower value. This is because our omega is
% smaller. 

%% RK4 Solver for the SIR Model
function [S, I, R] = RK4_SIR(beta_func, gamma, t, N)
    h = 0.1;
    n = length(t); % Number of time steps
    
    % Initialize S, I, R vectors
    S = zeros(n, 1);
    I = zeros(n, 1);
    R = zeros(n, 1);

    % Initial conditions
    S(1) = 990;
    I(1) = 10;
    R(1) = 0;

    % Solve the SIR equations using RK4
    for i = 1:(n-1)
        current_beta = beta_func(t(i)); % Evaluate beta at current time step

        k1_S = -current_beta * S(i) * I(i) / N;
        k1_I = current_beta * S(i) * I(i) / N - gamma * I(i);
        k1_R = gamma * I(i);

        k2_S = -beta_func(t(i) + h/2) * (S(i) + h * k1_S / 2) * (I(i) + h * k1_I / 2) / N;
        k2_I = beta_func(t(i) + h/2) * (S(i) + h * k1_S / 2) * (I(i) + h * k1_I / 2) / N - gamma * (I(i) + h * k1_I / 2);
        k2_R = gamma * (I(i) + h * k1_I / 2);

        k3_S = -beta_func(t(i) + h/2) * (S(i) + h * k2_S / 2) * (I(i) + h * k2_I / 2) / N;
        k3_I = beta_func(t(i) + h/2) * (S(i) + h * k2_S / 2) * (I(i) + h * k2_I / 2) / N - gamma * (I(i) + h * k2_I / 2);
        k3_R = gamma * (I(i) + h * k2_I / 2);

        k4_S = -beta_func(t(i) + h) * (S(i) + h * k3_S) * (I(i) + h * k3_I) / N;
        k4_I = beta_func(t(i) + h) * (S(i) + h * k3_S) * (I(i) + h * k3_I) / N - gamma * (I(i) + h * k3_I);
        k4_R = gamma * (I(i) + h * k3_I);

        % Update values
        S(i+1) = S(i) + (h/ 6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
        I(i+1) = I(i) + (h/ 6) * (k1_I + 2*k2_I + 2*k3_I + k4_I);
        R(i+1) = R(i) + (h/ 6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
    end
end
