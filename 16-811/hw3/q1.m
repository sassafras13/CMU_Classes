% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3
% Problem 1

%% 

clear all ; clc ; close all ; 

%% 

% sinh(0)
% 
% cosh(0)
% 
% factorial(5)

x = -3:0.1:3 ; 
fx = (1/3) + 2*sinh(x) ; 
y = -0.5 + 0.1.*x + 2.*x.^2 ; 

fig1 = figure(1)
plot(x, fx, '-b') ; 
hold on 
plot(x, y, '-g') ; 
title('Function f(x)') ; xlabel('x') ; ylabel('f(x)') ; 
