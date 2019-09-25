% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 6

%% 

clc ; clear all ; close all ; 

Q = [1, -4, 6, -4, 0 ; 
     0, 1, -4, 6, -4 ; 
     1, 2, -8, 0, 0 ; 
     0, 1, 2, -8, 0 ; 
     0, 0, 1, 2, -8 ] ; 
 
detQ = det(Q)

Ai = [-4 , 6, -4, 0 ; 
      1, -4, 6, -4 ; 
      2, -8, 0, 0 ; 
      1, 2, -8, 0 ] ; 
  
Aj = [1, 6, -4, 0 ; 
      0, -4, 6, -4 ; 
      1, -8, 0, 0 ; 
      0, 2, -8, 0 ] ; 
  
x = (-1)^3 * det(Ai) / det(Aj) 