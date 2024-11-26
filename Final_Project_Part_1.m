%% Final Project - Part 1

% Seasonal Influenza
beta1 = 0.3; % Transmission Rate
gamma1 = 0.1; % Recovery Rate

% COVID
beta2 = 1; % Transmission Rate
gamma2 = 0.1; % Recovery Rate

% Measles
beta3 = 2; % Transmission Rate
gamma3 = 0.2; % Recovery Rate



h = 1; % time step, days
N = 1000; % Total Population, constant
t = 0:h:100; % total simulation time, 100 days
S0 = 990; % Initial Condition for S(t)
I0 = 10; % Initial Condition for I(t)
R0 = 0; % Initial Condition for R(t)
y0 = [S0 I0 R0]; % initial conditions vector



% Model for Seasonal Inluenza, COVID, and Measles

SIR_influenza = @(t,y) [  -beta1*y(1) * y(2) / N; % dS 
                 beta1*y(1) * y(2) / N - gamma1 * y(2); % dI
                 gamma1 * y(2); % dR 
                                    ];

SIR_covid = @(t,y) [  -beta2*y(1) * y(2) / N; % dS 
                 beta2*y(1) * y(2) / N - gamma2 * y(2); % dI
                 gamma2 * y(2); % dR 
                                    ];

SIR_measles = @(t,y) [  -beta3*y(1) * y(2) / N; % dS 
                 beta3*y(1) * y(2) / N - gamma3 * y(2); % dI
                 gamma3 * y(2); % dR 
                                    ];


% S O L V E - Seasonal Influenza
y = RK4(SIR_influenza, t, y0);

% P L O T S - Seasonal Influenza
figure(1)
plot(t,y(:,1),'r-');
hold on
plot(t,y(:,2),'g-')
plot(t,y(:,3),'b-')
title('Influenza Simulation');
xlabel('time, t (days)');
ylabel('population');
legend('S(t)','I(t)','R(t)')



% S O L V E - COVID
y = RK4(SIR_covid, t, y0);

% P L O T S - COVID
figure(2)
plot(t,y(:,1),'r-');
hold on
plot(t,y(:,2),'g-')
plot(t,y(:,3),'b-')
title('COVID Simulation');
xlabel('time, t (days)');
ylabel('population');
legend('S(t)','I(t)','R(t)')



% S O L V E - Measles
y = RK4(SIR_measles, t, y0);

% P L O T S - Measles
figure(3)
plot(t,y(:,1),'r-');
hold on
plot(t,y(:,2),'g-')
plot(t,y(:,3),'b-')
title('Measles Simulation');
xlabel('time, t (days)');
ylabel('population');
legend('S(t)','I(t)','R(t)')



%% Runge-Kutta Method; 4th Order 
function y = RK4(f, t, y0)
h = 1;
y = zeros(length(t), length(y0)); % Initialized Solution Matrix
y(1,:) = y0;

for i = 1:(length(t)-1)
    k1 = f( t(i), y(i,:) )*h;
    k2 = f( t(i) + (h/2), y(i,:) + (k1/2) )*h;
    k3 = f( t(i) + (h/2), y(i,:) + (k2/2) )*h;
    k4 = f( t(i) + h, y(i,:) + k3 )*h;

    y(i+1,:) = y(i,:) + ( (h/6) * (k1' + 2*k2' + 2*k3' + k4') );

end
end
