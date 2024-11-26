
%% Final Project - Part 2
clc,  clear all 
% Seasonal Influenza
beta1 = 0.3; % Transmission Rate
gamma1 = 0.1; % Recovery Rate

N = 1000; % Total Population, constant
 % total simulation time, 100 days
S0 = 990; % Initial Condition for S(t)
I0 = 10; % Initial Condition for I(t)
R0 = 0; % Initial Condition for R(t)
y0 = [S0 I0 R0]; % initial conditions vector



% Model for Seasonal Inluenza, COVID, and Measles

SIR_influenza = @(t,y) [  -beta1*y(1) * y(2) / N; % dS 
                 beta1*y(1) * y(2) / N - gamma1 * y(2); % dI
                 gamma1 * y(2); % dR 
  ]


% S O L V E - Seasonal Influenza
y = RK4(SIR_influenza, y0, 2);
y1 = RK4(SIR_influenza, y0, 1); 





%% Interpolation ------------------------------------------------------------------------------------
n = 1; 
ODD  = 3:2:100; 
for  j = 1:3
    for n= 3:(51-2)
    x0  = n-1
    x1 = n+1
    x2 = n+2
    
        
        
    b0 = y((x0),j)
    b1 = (y((x1),j)-y((x0),j))./((x1)-(x0))
    b2 = ((y((x2),j)-y((x1),j))./((x2)-(x1)) - (y((x1),j)-y((x0),j))./((x1)-(x0)))./(x2-x0)
    
    LN(n,j) = b0 + b1*(n-x0) % Linear Newton
    QN(n,j) = b0 + b1*(n-x0)+b2*(n-x0)*(n-x1) %Quadratic Newton
    
    end
   E1(j) = sqrt(sum(LN(:,j) - y1(ODD,j)).^2) ./ (length(ODD)); %Error linear
    E2(j) = sqrt(sum(QN(:,j) - y1(ODD,j)).^2) ./ (length(ODD)); %Error Quadratic 
    
end

ERROR = [E1 ; E2]
array2table(ERROR,"VariableNames", {'S(t)', 'I(t)','R(t)'}, "RowNames", {'Linear','Quadratic'})


%% Runge-Kutta Method; 4th Order
function y = RK4(f,y0,h)
t = 0:h:100;
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



