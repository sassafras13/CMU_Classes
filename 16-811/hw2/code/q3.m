% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 3
% References: Dawkins, Paul. "Section 4-13: Newton's Method." Paul's Online Notes.
% <http://tutorial.math.lamar.edu/Classes/CalcI/NewtonsMethod.aspx> Visited
% 09/19/2019.

%% 
clc ; clear all ; close all ; 

%% Plot function to understand how it looks
% plot tan(x) - x = fx just to see what it looks like
x = 0:0.1:10 ; 
y = tan(x) - x ; 
line = zeros(length(x),1) ; 
seven = 7*ones(length(x),1) ; 
yseven = linspace(-80,80,length(x)) ; 

f = figure(1) ; 
plot(x,y,'-b') ; 
hold on
plot(x,line,'-r'); 
hold on 
plot(seven,yseven,'-g');
hold on
xlabel('x') ; ylabel('y') ; 

%% Bisection and Newton's method implementations

% write a function that takes in a starting point x0 and looks for the
% solution that converges to k decimal places. It should return xf.

% we want one solution less than 7 and one solution greater than 7 such
% that there are no other solutions in between 

% looking at the plot from above, I can see that the roots for this
% function that meet these requirements should be around x0_low = 4.5 and
% x0_high = 7.5

eps = 0.1 ; % accuracy required from bisection
k = 3 ; % number of decimal places of convergence for Newton's method

% reference point
xref = 7 ; 

% 2 starting values
x0_low = 6; 
x0_high = 7; 

% x0 = [x0_low, x0_high] ; 
x0 = x0_high ; 

xf = zeros(2,1) ; 

for i = 1:length(x0)
    % use bisection to get into the region of convergence
    xf_bisection = bisection(x0(i),eps) ;

    % use Newton's method to refine answer
    format long
    xf_newton = newtonmethod(xf_bisection,k) ; 
    
    xf(i) = xf_newton ; 
end

xf

figure(f)
plot(xf(1),0,'om','MarkerSize',10,'MarkerFaceColor','m') ; 
hold on 
plot(xf(2),0,'om','MarkerSize',10,'MarkerFaceColor','m') ; 
hold off

% check values are very high because the roots are right near poles of the
% function, indicates I have found the roots correctly
check_low = fx(xf(1)) 

check_high = fx(xf(2))

%% Functions

function [xf,i] = newtonmethod(x0,k)
    xn = 0 ; 
    xn1 = x0 ; 
    
    eps = 5*10^(-(k+1)) ; 
    
    i = 0 ; 
    
    while abs(xn1 - xn) > eps
        xn = xn1 ;
        xn1 = xn - (fx(xn) / fpx(xn)) ; 
        i = i + 1 ; % i = # iterations
        check = abs(xn1 - xn) ; 
    end
    
    xf = xn1 ; 
end

function y = fx(x)
    y = tan(x) - x ; 
end

function y = fpx(x)
    y = (sec(x))^2 - 1 ; 
end

function xf = bisection(x0, eps)

    % pick bracket [a,b] s.t. a > 0 and b < 0 or a < 0 and b > 0 
    % initial bracket
    xa = x0 - 5  
    xb = x0 + 5  
    ya = fx(xa)  
    yb = fx(xb)  
    
    % if the brackets are the same sign, move them until a > 0 and b < 0 
    if ya < 0 && yb < 0 
        while yb < 0
            xb = xb - eps ; 
            yb = fx(xb) ; 
        end
    elseif yb > 0 && ya > 0 
        while ya > 0
            xa = xa + eps ; 
            ya = fx(xa) ; 
        end
    end

    while abs(xa - xb) > eps
          
        % calculate midpoint (a + b)/2
        xmid = (xa + xb) / 2  
        ymid = fx(xmid)  
        
        % assign midpoint to a or b depending on sign
        if ymid > 0 && ya > 0
            xa = xmid ; 
        elseif ymid > 0 && yb > 0 
            xb = xmid ; 
        elseif ymid < 0 && ya < 0
            xa = xmid ; 
        elseif ymid < 0 && yb < 0
            xb = xmid ; 
        end
        
    end
    
    xf = (xa + xb) / 2 ; 
end
