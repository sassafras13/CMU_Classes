% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3, Resubmit 1
% Problem 1

%% 

clear all ; clc ; close all ; 

%% 

x = -3:0.1:3 ; 
fx = (1/3) + 2*sinh(x) ;
x0 = zeros(length(x)) ; 
interp = 5.3872*x + (1/3) ; 

fig1 = figure(1) ;
plot(x, fx, '-b') ; 
hold on 
plot(x, interp, '--r') ; 
title('Function f(x)') ; xlabel('x') ; ylabel('f(x)') ; 
legend('f(x)', 'p(x)') ; 

ex_inf = abs(2*sinh(x) - 5.3872*x) ; 
ex_inf = max(ex_inf) 


