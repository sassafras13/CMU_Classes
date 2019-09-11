
%% Problem 2
clc ; 

% matrix A1 from Problem 2
A1 = [10, 9, 2 ; 
     5, 3, 1 ; 
     2, 2, 2 ] ; 

% matrix A2 
A2 = [16, 16, 0, 0 ; 
     4, 0, -2, 0 ; 
     0, 1, -1, 0 ; 
     0, 0, 0, 1 ; 
     0, 0, 1, 1 ] ; 
 
% matrix A3 from Problem 2
A3 = [10, 6, 4 ; 
     5, 3, 2 ; 
     1, 1, 0 ] ; 

[U,S,V] = svd(A1) ;

%% Problem 3 Part (a)
clc ; 

b1 = [-2 ; 2 ; 4 ] ; 

[U1,S1,V1] = svd(A1) ;

invS1 = zeros(size(S1,1)) ; 

for i = 1:size(S1,1)
    invS1(i,i) = 1/S1(i,i) ; 
end

invS1 ;

invA1 = inv(A1)  % it is invertible and square
nullA1 = null(A1) % it has a trivial null space

% A1 has an exact solution, xbar_A1
xbar_A1 = V1*invS1*U1'*b1 

%% Problem 3 Part (b)
clc ; 

b2 = [2 ; 1; -1 ] ; 

[U2,S2,V2] = svd(A3) 

invS2 = zeros(size(S2,1)) ; 
eps = 10^-10 ; 

for i = 1:size(S2,1)
    if S2(i,i) > eps
        invS2(i,i) = 1/S2(i,i) ; 
    else 
        invS2(i,i) = 0 ; 
    end
end

invA3 = inv(A3) % A3 is non-invertible 
nullA3 = null(A3) % A3 has a non-trivial null space

% A3 has multiple solutions of the form xbar_A3 + xN
xbar_A3 = V2*invS2*U2'*b2 

%% Problem 3 Part (c)
clc ; 

b3 = [1 ; 3 ; -1 ] ; 

% A3 has multiple solutions of the form xbar_A3_b3 + xN
xbar_A3_b3 = V2*invS2*U2'*b3

%% Problem 5
clc ;

% let A1 represent 3 points in 3D space (each row is a point)

% let's apply a rotation about each axis
% syms theta
theta = pi/6 ; 

Rx = [1, 0, 0 ; 
      0, cos(theta), -sin(theta) ; 
      0, sin(theta), cos(theta) ] ; 
  
Ry = [cos(theta), 0, sin(theta) ; 
      0, 1, 0 ; 
      -sin(theta), 0, cos(theta) ] ; 
  
Rz = [cos(theta), -sin(theta), 0 ; 
      sin(theta), cos(theta), 0 ; 
      0, 0, 1 ] ; 
  
R = Rz*Ry*Rx ; 

% and let's apply a translation
% T = [1, 0, 0, 2 ; 
%      0, 1, 0, 3 ; 
%      0, 0, 1, -2 ; 
%      0, 0, 0, 1 ] ; 

v = [2 ; 3 ; -2 ] ; 

A = [R , v ; zeros(1,3), 1]  

x = [A1, zeros(3,1) ; zeros(1,3) , 1 ] ; 

b = A*x  

checkA = b*inv(x) 