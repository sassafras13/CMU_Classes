% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 7
% References: 

%% 
clc ; clear all ; close all ; 

%% Plot function to understand how it looks

x = -5:0.1:10 ; 
y = -5:0.1:10 ; 
% zp = px(x,y) ;  
% zq = qx(x,y) ;

yp = pxZeroContourY(x) ; 
yq = qxZeroContourY(x) ; 
xp = pxZeroContourX(y) ; 
xq = qxZeroContourX(y) ; 

f = figure(1) ; 
plot(xp,yp) ; 
hold on
plot(xq,yq) ;  
hold on 
xlabel('x') ; ylabel('y')  
grid on 
legend('p','q') ; 

%% Symbolically solve for the determinant of Q

syms x 'real'

Q = [2 , 4, (3 + 2*(x^2) + 4*x), 0 ; 
     0, 2, 4, (3 + 2*(x^2) + 4*x) ; 
     1, (2*x + 5), (4 + (x^2) + 3*x), 0 ; 
     0, 1, (2*x + 5), (4 + (x^2) + 3*x) ]  
 
resultant = det(Q)

C = [16,80,144,100,19] ; 
roots(C)

%% Functions 

function zp = px(x,y) 
    zp = 2.*(x.^2) + 2.*(y.^2) + 4.*x + 4.*y + 3 ; 
end

function zq = qx(x,y)
    zq = x.^2 + y.^2 + 2.*x.*y + 3.*x + 5.*y + 4 ; 
end

function yp = pxZeroContourY(x)
    
    yp = zeros(length(x),1) ; 
    
    for i = 1:length(x)
        C = [2,4,(3 + 2*(x(i)^2) + 4*x(i))] ;
        rootC = roots(C) ;
        yp(i) = real(rootC(1)) ;
    end
   
end

function yq = qxZeroContourY(x)
    
    yq = zeros(length(x),1) ; 
    
    for i = 1:length(x)
        C = [1,(2*x(i) + 5),(4 + (x(i)^2) + 3*x(i))] ; 
        rootC = roots(C) ;
        yq(i) = real(rootC(1)) ;
    end
   
end

function xp = pxZeroContourX(y)
    
    xp = zeros(length(y),1) ; 
    
    for i = 1:length(y)
        C = [2,4,(3 + 2*(y(i)^2) + 4*y(i))] ; 
        rootC = roots(C) ;
        xp(i) = real(rootC(1)) ;
    end
   
end

function xq = qxZeroContourX(y)
    
    xq = zeros(length(y),1) ; 
    
    for i = 1:length(y)
        C = [1, (2*y(i) + 3), (4 + (y(i)^2) + 5*y(i))] ; 
        rootC = roots(C) ;
        xq(i) = real(rootC(1)) ;
    end
   
end