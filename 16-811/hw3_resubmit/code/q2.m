% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3, Resubmit 1
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

figure(1)
% plot(xi, fi,'-b') ; 
% hold on

% The plot shows that the function could really be divided into 2, at the
% point (x,y) = (3,27) which is at index i = 30

fna = fi(1:30) ; 
fnb = fi(30:end) ; 
xia = 1:1:30 ; 
xia = xia/10 ; 
xib = 30:1:100 ; 
xib = xib/10 ; 

% first function a
phia1 = ones(length(fna), 1) ; 
phia2 = xia' ; 
phia3 = (phia2).^2 ;
phia4 = (phia2).^3 ; 
phia5 = (phia2).^4 ; 

% second function b
phib1 = ones(length(fnb), 1) ; 
phib2 = xib' ; 
phib3 = (phib2).^2 ;
phib4 = (phib2).^3 ; 
phib5 = (phib2).^4 ; 

Aa = [phia1, phia2, phia3, phia4] ; 

Ab = [phib1, phib2, phib3, phib4] ; 

% check matrix
% invA = inv(A) 
% nullA = null(A) 

[Ua,Sa,Va] = svd(Aa) ;
[Ub,Sb,Vb] = svd(Ab) ; 

% calculate inverse element by element
invSa = zeros(size(Sa,1),size(Sa,2)) ; 

for i = 1:size(Sa,1)
    for j = 1:size(Sa,2) 
            if Sa(i,j) > eps
            invSa(i,j) = 1/Sa(i,j) ; 
        else 
            invSa(i,j) = 0 ; 
            end
    end
end

xbara = Va*invSa'*Ua'*fna ;

invSb = zeros(size(Sb,1),size(Sb,2)) ; 

for i = 1:size(Sb,1)
    for j = 1:size(Sb,2) 
            if Sb(i,j) > eps
            invSb(i,j) = 1/Sb(i,j) ; 
        else 
            invSb(i,j) = 0 ; 
            end
    end
end

xbarb = Vb*invSb'*Ub'*fnb ;

pxa = xbara(1) + xbara(2)*xia + xbara(3)*(xia.^2) + xbara(4)*(xia.^3) ;
plot(xia, pxa, '-r') ;
hold on

pxb = xbarb(1) + xbarb(2)*xib + xbarb(3)*(xib.^2) + xbarb(4)*(xib.^3) ;
plot(xib, pxb, '-g') ; 

legend('Function A','Function B') ; title('Approximating Functions') ; 
