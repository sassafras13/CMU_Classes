% 16-811
% hw5 q3

%% 

clc ; clear all ; close all ; 

%% 

syms g yp y0 y gy gx 'real'

eqn = (1/sqrt(-2*g))*(yp/sqrt((y0-y)*(1+yp^2))) - (gy/(gx + gy*yp))*(1/(sqrt(-2*g)))*(sqrt(1 + yp^2)/sqrt(y0-y)) == 0 

Sa = solve(eqn,yp)

eqn2 = (yp^3)*(2*gx*gy) + (yp^2)*(gx^2 - 2*gy^2) - (gy^2) == 0

Sa2 = solve(eqn2, yp)

eqn3 = eqn2 / (yp - (gy/gx)) 

Sa3 = solve(eqn3, yp)