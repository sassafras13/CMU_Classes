% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3
% Problem 2

%% 
clear all ; close all ; clc ; 

%%
fid = fopen('problem2.txt') ; 
data = textscan(fid, '%f') ; 
fi = data{1,1} ; 
fi = fi(2:end) ; 

i = 1:1:100 ; 
xi = i/10 ; 

plot(xi, fi,'-b') ; 
hold on

phi1 = ones(length(fi), 1) ; 
phi2 = xi' ; 
phi3 = (phi2).^2 ;
phi4 = (phi2).^3 ; 

A = [phi1, phi2, phi3, phi4] ; 

% check matrix
% invA = inv(A) 
% nullA = null(A) 

[U,S,V] = svd(A) ;

% calculate inverse element by element
invS = zeros(size(S,1),size(S,2)) ; 

for i = 1:size(S,1)
    for j = 1:size(S,2) 
            if S(i,j) > eps
            invS(i,j) = 1/S(i,j) ; 
        else 
            invS(i,j) = 0 ; 
            end
    end
end

xbar = V*invS'*U*fi

px = -3.2033 + 20.9635*xi - 5.6690*(xi.^2) + 0.4086*(xi.^3) ; 
plot(xi, px, '-r') ; 

