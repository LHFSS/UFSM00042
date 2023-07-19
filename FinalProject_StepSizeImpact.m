clc; close all; clear;

global K;
K = 250;

global r;
r = 0.2;

N0 = 10;
h = 3.5;
h_A = 0.5;
ti = 0;
tf = 35;

for h = 0.1:0.1:5
    [t_R, N_R] = RungeKutta(ti, tf, h, N0);
    [t_E, N_E] = Euler(ti, tf, h, N0);
    [t_A, N_A] = EqLogAnalitica(ti, tf, h_A, N0);
    clf;
    plot(t_A, N_A);
    hold on;
    plot(t_R, N_R);
    hold on;
    plot(t_E, N_E);
    grid on;
    legend('Analítico', 'Runge-Kutta', 'Euler');
    pause(0.05);
end


function [t_R, N_R] = RungeKutta(ti, tf, h, N0)
    t_R = ti:h:tf+h;
    n = length(t_R);
    N_R = zeros([1 n]);
    N_R(1) = N0;
    
    for i = 1:n-1
        k1 = h * EqLog(t_R(i), N_R(i));
        k2 = h * EqLog(t_R(i) + h/2, N_R(i) + h*k1/2);
        k3 = h * EqLog(t_R(i) + h/2, N_R(i) + h*k2/2);
        k4 = h * EqLog(t_R(i) + h, N_R(i) + h*k3);
        
        N_R(i+1) = N_R(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end

function [t_E, N_E] = Euler(ti, tf, h, N0)
    t_E = ti:h:tf+h;
    n = length(t_E);
    N_E = zeros([1 n]);
    N_E(1) = N0;
    
    for i = 1:n-1
        N_E(i+1) = N_E(i) +  h * EqLog(t_E(i), N_E(i));
    end
end

function EqLog = EqLog(t, N)
    global r; global K;
    EqLog = r*N*((1 - N/K));
end

function [t_A, N_A] = EqLogAnalitica(ti, tf, h, N0)
    t_A = ti:h:tf+h;
    global r; global K;
    
    C = (1/N0) - (1/K);
    N_A = K ./ (1 + C*K*exp(-r.*t_A));
end
