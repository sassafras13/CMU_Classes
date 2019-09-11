% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 1
% Problem 3

%% 
clc ; clear all ; close all ; 

eps = 10^-10 ; % epsilon for checking near-zero values

%% Part (a)
A1 = [10 , 9 , 2 ; 
      5 , 3 , 1 ; 
      2 , 2 , 2 ] ; 
  
b1 = [-2 ; 2 ; 4 ] ; 

% check matrix
invA1 = inv(A1)  % it is invertible and square
nullA1 = null(A1) % it has a trivial null space

% given that matrix is invertible and has a trivial null space, we expect to find one exact
% solution using SVD decomposition

% calculate SVD solution 
[U1,S1,V1] = svd(A1) ;

% calculate inverse element by element
invS1 = zeros(size(S1,1)) ; 

for i = 1:size(S1,1)
    
    if S1(i,i) > eps
        invS1(i,i) = 1/S1(i,i) ; 
    else 
        invS1(i,i) = 0 ; 
    end
end

% calculate exact solution xbar 
xbar_A1 = V1*invS1*U1'*b1 

%% Part (b)
clc ; 

A2 = [10 , 6 , 4 ; 
      5 , 3 , 2 ; 
      1 , 1 , 0 ] ; 
  
b2 = [2 ; 1; -1 ] ; 

% check matrix
invA2 = inv(A2) % A2 is non-invertible 
nullA2 = null(A2) % A2 has a non-trivial null space

% given that matrix is non-invertible and has a non-trivial null space, we expect to find many exact
% solutions using SVD decomposition

% check column space of A2
symA2 = sym(A2) ; 
colspaceA2 = colspace(symA2) 

[U2,S2,V2] = svd(A2) 

% calculate inverse element by element
invS2 = zeros(size(S2,1)) ; 

for i = 1:size(S2,1)
    if S2(i,i) > eps
        invS2(i,i) = 1/S2(i,i) ; 
    else 
        invS2(i,i) = 0 ; 
    end
end

% A2 has one unique solution
xbar_A2 = V2*invS2*U2'*b2 

%% Problem 3 Part (c)
clc ; 

b3 = [1 ; 3 ; -1 ] ; 

% vector b3 is not in the column space of matrix A2 so the best solution
% we can find is the least squares solution (i.e. the SVD solution)
xbar_A2_b3 = V2*invS2*U2'*b3
