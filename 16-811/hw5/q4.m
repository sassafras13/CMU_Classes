% hw5 q4

syms theta1 theta2 l1 l2 g m1 m2 thetadot1 thetadot2 'real'

h1 = l1*sin(theta1) 
h2 = l1*sin(theta1) - l2*cos(theta1 + theta2) 
h3 = l1*sin(theta1) + l2*sin(theta1 + theta2) 

xdot = -l1*sin(theta1)*thetadot1 - l2*sin(theta1 + theta2)*(thetadot1 + thetadot2) 
ydot = l1*cos(theta1)*thetadot1 + l2*cos(theta1 + theta2)*(thetadot1 + thetadot2) 

v1 = l1*thetadot1 
v22 = xdot^2 + ydot^2 % v2 squared

