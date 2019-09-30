% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 1
% Problem 1

%% 
clc ; clear all ; close all ; 

% test matrices; uncomment to run code on desired matrix A
% matrix A1 from Problem 2
A = [10, 9, 2 ; 
     5, 3, 1 ; 
     2, 2, 2 ] ; 

% matrix A3 from Problem 2
% A = [10, 6, 4 ; 
%      5, 3, 2 ; 
%      1, 1, 0 ] ; 

% sample
% A = [2, 3, -1, 1 ; 
%      -4, -4, 3, -2 ; 
%      2, 9, 1, 4 ; 
%      4, 6, -1, 3 ] ; 
 
[L,D,U,P] = LUdecomp(A) 

PA = P*A

LDU = L*D*U 

%% Function definition

function [L,D,U,P] = LUdecomp(A) 
    %LUdecomp performs the LU decomposition on the matrix A. It assumes A
    %is square and invertible. This function was written by Emma
    %Benjaminson and draws on information provided by: 
    %https://www.csun.edu/~panferov/math262/262_rref.pdf
    %MATLAB function rref()
    
    [m, ~] = size(A) ; % assume A is square so just need one dimension
    
    i = 1 ; % row index 
    j = 1 ; % col index
    
    eps = 0.001 ; % some tolerance for small numbers that should be 0
    
    L = zeros(m,m) ; % matrix L is going to track Gaussian elimination process
    P = eye(m,m) ; % matrix P tracks permutations
        
    % first run a loop to zero out very small values and order rows in
    % decreasing size
    while i <= m && j <= m 
        
        % look for max non-zero pivot element in current col
        [M, k] = max(abs(A(i:m, j))) ; 
        k = k + i -1 ; % because k will be an index with respect to A(i:m, j), not A(:,:)
        
        % check if any value is small enough that it should be zero
        if M < eps
            A(i:m,j) = 0 ; 
            j = j + 1 ; 
        else 
            % swap current row with the row containing the max non-zero pivot
            % element
            A([i k], j:m) = A([k i], j:m) ; 
            
            % also perform the swap for L and P 
            L([i k], :) = L([k i], :) ; 
            P([i k], :) = P([k i], :) ; 
            
            for k = [i+1:m] % work through all subsequent rows
                            
                % divide the pivot by the element in the row above it (same
                % column) to get the factor x
                x = A(k, j)/A(i,j) ; 
            
                % use x to subtract the row above 
                A(k, j:m) = A(k, j:m) - x*A(i,j:m) ; 
            
                L(k, i) = x ;                 
            end
        end

        i = i + 1 ; 
        j = j + 1 ; 
        
    end
    
   Aprime = A ; 
    
    for i = 1:m
        L(i,i) = 1 ; 
        D(i,i) = Aprime(i,i) ; 
        U(i,i) = 1 ; 
        for j = i+1:m
            U(i,j) = Aprime(i,j) / D(i,i) ; 
        end
    end
end
