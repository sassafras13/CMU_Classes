% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3
% Problem 3

%% 
clear all ; close all ; clc ; 

%% 
x = -1:0.1:1 ; 

t3 = 4*x.^3 - 3*x ; 
t4 = 8*x.^4 - 8*x.^2 + 1 ; 

figure(1)
plot(x, t3, '-b') ; 
hold on 
plot(x, t4, '-g') ; 
xlabel('x') ; ylabel('y') ; title('Problem 3') ; 
legend('T3','T4') ; 