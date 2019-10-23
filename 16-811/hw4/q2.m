% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 4
% Problem 2
% References: 
%% 
clc ; clear all ; close all ; 

%% Part a

x = -2:0.01:2 ; 
y = -2:0.01:2 ; 
dfdx = fx(x) ; 
dfdy = fy(y) ; 

figure(1)
plot(x,dfdx,'-r') ; 
hold on
plot(y,dfdy,'-b') ; 
hold off
title('dfdx and dfdy') ; xlabel('x') ; ylabel('y') ; 
legend('dfdx','dfdy') ; 

%% Functions

function dfdx = fx(x)
    dfdx = 3*x.^2 - 4*x ; 
end

function dfdy = fy(y) 
    dfdy = 3*y.^2 + 6*y ; 
end