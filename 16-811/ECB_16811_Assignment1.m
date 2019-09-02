% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 1

%% 
clc ; clear all ; close all ; 

%% Problem 1

A = [1, -2, 1 ; 
     1, 2, 2 ; 
     2, 3, 4 ]  

% A = [1, 1, 0 ; 
%      1, 1, 2 ; 
%      4, 2, 3 ]  
 
[L,D,U,Aprime] = LUdecomp(A) 