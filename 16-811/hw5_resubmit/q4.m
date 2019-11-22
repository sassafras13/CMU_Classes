% 16-811
% hw5 q4

%% 
clear all ; close all ; clc ; 

syms l1 l2 t1 t2 tdot1 tdot2 'real'

xdot2 = -l1*sin(t1)*tdot1 - l2*sin(t1 + t2)*(tdot1 + tdot2) ; 
ydot2 = l1*cos(t1)*tdot1 + l2*cos(t1 + t2)*(tdot1 + tdot2) ; 

v22 = (xdot2)^2 + (ydot2)^2 

simplify(v22)

xdot3 = -l1*sin(t1)*tdot1 + l2*sin(t1 + t2)*(tdot1 + tdot2) ; 
ydot3 = l1*cos(t1)*tdot1 - l2*cos(t1 + t2)*(tdot1 + tdot2) ; 

v32 = simplify((xdot3)^2 + (ydot3)^2) 