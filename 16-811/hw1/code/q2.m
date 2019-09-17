% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 1
% Problem 2

%% 
clc ; clear all ; close all ; 

% test matrices; uncomment to run code on desired matrix A
% matrix A1
% A = [10, 9, 2 ; 
%      5, 3, 1 ; 
%      2, 2, 2 ] ; 

% matrix A2 
% A = [16, 16, 0, 0 ; 
%      4, 0, -2, 0 ; 
%      0, 1, -1, 0 ; 
%      0, 0, 0, 1 ; 
%      0, 0, 1, 1 ] ; 
 
% matrix A3 from Problem 2
A = [10, 6, 4 ; 
     5, 3, 2 ; 
     1, 1, 0 ] ; 

[U,S,V] = svd(A) 