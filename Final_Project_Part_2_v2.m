%% Final Project - Part 1
clc, clear all 
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

%% Interpolation ------------------------------------------------------------------------------------

%Initializing
ODD  = 3:2:100; 
S = S_influenza
I = I_influenza
R = R_influenza
j= 0



for j = 1:length(S)-1 %running till I have all values
    j=j+1 %Counter
    
    
    if (-1)^j > 0       % (-1)^ even = 1 therefore these values will stay the same
        S_Linear_Newton(j) = S(j)
        S_Quad_Newton(j) = S(j)
        I_Linear_Newton(j) = I(j)
        I_Quad_Newton(j) = I(j)
        R_Linear_Newton(j) = R(j)
        R_Quad_Newton(j) = R(j)
        
        
    end
    
    if (-1)^j < 0   % (-1)^Odd = -1 therefore these values will be interpolated
        if j < 3 % Can't do x0 = n-1 when n = 1 therefore substitute these values
            x0  = j+1
            x1 = j+2
            x2 = j+3
        elseif j>97 % Can't do x2 = n+2 when n = 101 because there's no value of SIR(101+2)
            x0 = j-1
            x1 = j-2
            x2 = j-3
        else %General Newton form
            x0  = j-1
            x1 = j+1
            x2 = j+2
        end
        
        [S_Linear_Newton(j) , S_Quad_Newton(j)]= Interpolation(S,x0,x1,x2,j)
        [I_Linear_Newton(j) , I_Quad_Newton(j)]= Interpolation(I,x0,x1,x2,j)
        [R_Linear_Newton(j) , R_Quad_Newton(j)]= Interpolation(R,x0,x1,x2,j)
        
        
    end
    
        
        
        
end   


Error  = @(f,F2) sqrt( sum( ( f'-F2).^2)  )/ length(F2)

S_Linear_Error = Error(S_Linear_Newton,S)
S_Quad_Error = Error(S_Quad_Newton,S)    
I_Linear_Error = Error(I_Linear_Newton,I)    
I_Quad_Error = Error(I_Quad_Newton,I)
R_Linear_Error = Error(R_Linear_Newton,R)
R_Quad_Error = Error(R_Quad_Newton,R)
 
Linear_ERROR = [S_Linear_Error, I_Linear_Error , R_Linear_Error]
Quad_ERROR = [S_Quad_Error, I_Quad_Error , R_Quad_Error]

ERROR = [Linear_ERROR ; Quad_ERROR]
array2table(ERROR,"VariableNames", {'S(t)', 'I(t)','R(t)'}, "RowNames", {'Linear','Quadratic'})

%% Interpolation Solver
function [Linear, Quadratic] = Interpolation(y,x0,x1,x2,j)
    b0 = y((x0))
    b1 = (y(x1)-y(x0))./((x1)-(x0))
    b2 = ((y(x2)-y(x1))./((x2)-(x1)) - (y(x1)-y(x0))./((x1)-(x0)))./(x2-x0)
       
   Linear = b0 + b1*(j-x0) % Linear Newton
    Quadratic = b0 + b1*(j-x0)+b2*(j-x0)*(j-x1) %Quadratic Newton
end




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
