% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 4 - Resubmission 1
% Problem 1
% References: Drakos, Nikos, Moore, Ross, Robert. "Fourth Order Runge Kutta
% Method" 2002-01-28.
% <http://www.math.ubc.ca/~israel/m215/runge/runge.html> Visited
% 10/20/2019.

%% 
clc ; clear all ; close all ; 

%% Part b - Euler's Method

% defining the diff eq 
h = 0.05 ; % step size
xi = 2:-h:1 ; % interval for x
y0 = sqrt(2) ; % initial value

% Euler's method
yfx = fx(xi) ; % true values

yi_euler = euler(xi,y0,h) ; % Euler's method estimate

e = abs(yfx - yi_euler') ; 

figure(1) 
subplot(2,1,1) 
plot(xi, yfx, '-or') 
hold on 
plot(xi, yi_euler, '-ob') 
axis([1 2 0 1.5]) ; 
xlabel('X') ; ylabel('Y') ; title('Comparing Eulers Method to True Value') ; 
legend('True Values','Eulers Method') ; 

figure(1) 
subplot(2,1,2)
plot(xi, e, '-ob') 
axis([1 2 0 0.5]) ; 
xlabel('X') ; ylabel('Error') ; title('Error in Eulers Method over [1,2]') ; 

[xi', yfx', yi_euler, e']

%% Part c - Runge-Kutta

yi_rk = rungekutta(xi,y0,h) ; 

e = abs(yfx - yi_rk') ; 

figure(2) 
subplot(2,1,1) 
plot(xi, yfx, '-or') 
hold on 
plot(xi, yi_rk, '-ob') 
axis([1 2 0 1.5]) ; 
xlabel('X') ; ylabel('Y') ; title('Comparing 4th Order Runge-Kutta to True Value') ; 
legend('True Values','4th Order Runge-Kutta','Location','southeast') ; 

figure(2)
subplot(2,1,2)
plot(xi, e, '-ob') 
axis([1 2 0 0.5]) ;
xlabel('X') ; ylabel('Error') ; title('Error in 4th Order Runge-Kutta over [1,2]') ; 

[xi', yfx', yi_rk, e']

%% Part d - Adams-Bashforth

yi_ab = adamsbashforth(xi,h) ; 

e = abs(yfx - yi_ab(4:end)') ; 

figure(3) 
subplot(2,1,1)
plot(xi, yfx, '-or') 
hold on 
plot(xi, yi_ab(4:end), '-ob') 
axis([1 2 0 1.5]) ; 
xlabel('X') ; ylabel('Y') ; title('Comparing 4th Order Adams-Bashforth to True Value') ; 
legend('True Values','4th Order Adams-Bashforth') ; 

figure(3)
subplot(2,1,2)
plot(xi, e, '-ob') 
axis([1 2 0 0.5]) ;
xlabel('X') ; ylabel('Error') ; title('Error in 4th Order Adams-Bashforth over [1,2]') ; 

[xi', yfx', yi_ab(4:end), e']

%% FUNCTIONS 

% true solution
function yfx = fx(xi) 
    yfx = sqrt(2) * sqrt(xi - 1) ; 
end

% Euler's method
function yi = euler(xi,y0,h) 
    yi = zeros(length(xi),1) ; 
    yi(1) = y0 ; 
    
    for i = 1:(length(xi)-1)
        yi(i+1) = yi(i) - h*(1/yi(i)) ; 
    end
end

% Runge-Kutta 4th order
function yi = rungekutta(xi,y0,h) 
    yi = zeros(length(xi),1) ; 
    yi(1) = y0 ; 
    
    for i = 1:(length(xi)-1)
        k1 = h*(1/yi(i)) ; 
        k2 = h*((1/yi(i)) + (k1/2)) ; 
        k3 = h*((1/yi(i)) + (k2/2)) ;
        k4 = h*((1/yi(i)) + k3) ; 
        yi(i+1) = yi(i) - (1/6)*(k1 + 2*k2 + 2*k3 + k4) ; 
    end
end

% Adams-Bashforth 4th order
function yi = adamsbashforth(xi,h)
    yi = zeros(length(xi)+3,1) ; 
    
    % starting values provided in question
    yi(1) = 1.51657508881031 ; 
    yi(2) = 1.48323969741913 ; 
    yi(3) = 1.44913767461894 ; 
    yi(4) = 1.4142135623731 ; 
    
    for i = 4:(length(yi) -1)
        fn3 = 1/yi(i-3) ; 
        fn2 = 1/yi(i-2) ; 
        fn1 = 1/yi(i-1) ; 
        fn = 1/yi(i) ; 
        
        yi(i+1) = yi(i) - (h/24)*(55*fn - 59*fn1 + 37*fn2 - 9*fn3) ; 
    end
    
end