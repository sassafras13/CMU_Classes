% 24-352 Spring 2020
% HW2
% Emma Benjaminson

%% 

clear all ; close all ; clc ; 

%% Q1

t = (0:0.01:2.5)' ; 
y = (t.^2) ./ 200 ;

s = tf('s') ; 
Ys = 1 / (100 * (s^2)) ; 
[Yt, tstep] = step(Ys,t) ; 

figure(1)
plot(t, y, '-b') 
hold on
plot(tstep, Yt, '--r') 
legend('Analytic Solution','Step Fn') 

%% Q2 

l1 = 0.75 ; 
l2 = 1.65 ; 
k = 420 ; 
theta0 = 0.03 ; 
omega0 = 0.015 ; 
m = 2.5 ; 
g = 9.8 ; 

omegan = sqrt( ( (m*g*(l1 + l2)) + k*( (l1^2) + (l1 + l2)^2 ) ) / (m * (l1 + l2)^2) ) ; 
T = (2*pi) / omegan ; 

t = 0:0.01:(5*T) ; 
y = theta0 * cos(omegan*t) + (omega0 / omegan) * sin(omegan*t) ; 

figure(2) 
plot(t,y,'-b')




