% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 4
% Problem 1
% References: 

%% 
clc ; clear all ; close all ; 

%% Part b - Euler's Method

% defining the diff eq 
xi = 2:-0.05:1 ; % interval for x
y0 = sqrt(2) ; % initial value
h = 0.05 ; % step size

% Euler's method
yfx = fx(xi) ; % true values

yi_euler = euler(xi,y0,h) ; % Euler's method estimate

e = abs(yfx - yi_euler') ; 

figure(1) 
plot(xi, yfx, '-or') 
hold on 
plot(xi, yi_euler, '-ob') 
xlabel('X') ; ylabel('Y') ; title('Comparing Eulers Method to True Value') ; 
legend('True Values','Eulers Method') ; 

figure(2) 
plot(xi, e, '-ob') 
xlabel('X') ; ylabel('Error') ; title('Error in Eulers Method over [1,2]') ; 

%% Part c - Runge-Kutta

yi_rk = rungekutta(xi,y0,h) ; 

e = abs(yfx - yi_rk') ; 

figure(3) 
plot(xi, yfx, '-or') 
hold on 
plot(xi, yi_rk, '-ob') 
xlabel('X') ; ylabel('Y') ; title('Comparing 4th Order Runge-Kutta to True Value') ; 
legend('True Values','4th Order Runge-Kutta') ; 

figure(4) 
plot(xi, e, '-ob') 
xlabel('X') ; ylabel('Error') ; title('Error in 4th Order Runge-Kutta over [1,2]') ; 

%% Part d - Adams-Bashforth

%% FUNCTIONS 

% true solution
function yfx = fx(xi) 
    yfx = log(xi) + 0.7211 ; 
end

% Euler's method
function yi = euler(xi,y0,h) 
    yi = zeros(length(xi),1) ; 
    yi(1) = y0 ; 
    
    yn = y0 ; 
    for i = 2:length(xi)
        yi(i) = yn - h*(1/xi(i)) ; 
        yn = yi(i) ; 
    end
end

% Runge-Kutta 4th order
function yi = rungekutta(xi,y0,h) 
    yi = zeros(length(xi),1) ; 
    yi(1) = y0 ; 
    
    yn = y0 ; 
    for i = 2:length(xi)
        k1 = h*yn ; 
        k2 = h*((1/(xi(i) + (h/2))) + (k1/2)) ; 
        k3 = h*((1/(xi(i) + (h/2))) + (k2/2)) ;
        k4 = h*((1/(xi(i) + h)) + k3) ; 
        yi(i) = yn - (1/6)*(k1 + 2*k2 + 2*k3 + k4) ; 
        yn = yi(i) ; 
    end
end