% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3
% Problem 1

%% 

clear all ; clc ; close all ; 

%% 

x = -3:0.1:3 ; 
fx = (1/3) + 2*sinh(x) ;
fx2 = 2*sinh(x) + 20 ; 
x0 = zeros(length(x)) ; 
interp = 3.0862*x + (1/3) ; 

fig1 = figure(1)
plot(x, fx, '-b') ; 
hold on 
plot(x, interp, '--r') ; 
title('Function f(x)') ; xlabel('x') ; ylabel('f(x)') ; 
legend('f(x)', 'p(x)') ; 

ex0 = 2*sin(-3) + (1/3) - (3.0862*-3 + 1/3)
ex1 = 2*sin(-1) + (1/3) - (3.0862*-1 + 1/3)
ex2 = 2*sin(1) + (1/3) - (3.0862*1 + 1/3)
ex3 = 2*sin(3) + (1/3) - (3.0862*3 + 1/3)

ex_inf = abs(2*sinh(x) - 3.0862*x) ; 
ex_inf = max(ex_inf) 


