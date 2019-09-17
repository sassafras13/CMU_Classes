% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 1
% Problem 5

%% 
clc ; clear all ; close all ; 

%% definitions

% matrix A is a 3 x 3 matrix where each column is a point, so each row is
% [x, y, z] coordinates for that point

%% create test points

% let A represent 3 points in 3D space (each row is a point) before
% transformation
A = rand(3,3) ; 

% let's apply a rotation about each axis
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
  
Rtest = Rz*Ry*Rx ; 

% let's produce a random translation
ttest = rand(3,1) ; 

B = Rtest*A + ttest ; 

% let's plot the points before and after
figure(1) 
scatter3(A(1,:),A(2,:),A(3,:)) ; 
labelsA = {'1','2','3'} ; 
text(A(1,:),A(2,:),A(3,:),labelsA) ; 
hold on 
scatter3(B(1,:),B(2,:),B(3,:)) ; 
labelsB = {'1"','2"','3"'} ; 
text(B(1,:),B(2,:),B(3,:),labelsB) ; 
hold off 
legend('A (before transform)','B (after transform)') ; 

%% extract rotation and translation matrices

[R, t] = transformation(A,B) ; 

% print and compare matrices
R 

Rtest

t

ttest

%% transformation

function [R, t] = transformation(A, B) 

    %transformation assumes that the input matrices A and B are square

    %%%%%
    % A %
    %%%%%
    
    % calculate centroid
    nA = size(A,1) ;
    Acen = (1/nA) * sum(A,2) ; % take the average of all the rows for the average centroid
    
    % find the vectors v that are the distance from each point to the
    % centroid
    vA = zeros(nA,nA) ; 
    
    for i = 1:nA
        vA(i,:) = A(:,i) - Acen ;
    end
    
    %%%%%
    % B %
    %%%%%
    
    % calculate centroid
    nB = size(B,1) ;
    Bcen = (1/nB) * sum(B,2) ; % take the average of all the cols for the average centroid
    
    % find the vectors v that are the distance from each point to the
    % centroid
    vB = zeros(nB,nB) ; 
    
    for i = 1:nB
        vB(i,:) = B(:,i) - Bcen ;
    end 
    
    %%%%%%%%%%%
    % H, R, t %
    %%%%%%%%%%%
    
    H = vA*vB' ;
    
    [U,~,V] = svd(H) ; 
    
    X = V*U' ; 
    
    R = X ; 
    
    t = B(:,1) - R*A(:,1) ;    
%     t2 = B(:,2) - R*A(:,2)  
%     t3 = B(:,3) - R*A(:,3)  
end
