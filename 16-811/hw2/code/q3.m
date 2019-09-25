% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 3
% References: Dawkins, P. "Section 4-13: Newton's Method." Paul's Online Notes.
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

eps = 0.001 ; % accuracy required from bisection
k = 3 ; % number of decimal places of convergence for Newton's method

% reference point
xref = 7 ; 

[xf_low,xf_high] = bisection(xref,eps) ;

% use Newton's method to refine answer
format long
xf_low_newton = newtonmethod(xf_low,k) ; 
xf_high_newton = newtonmethod(xf_high,k) ; 

figure(f)
plot(xf_low_newton,0,'om','MarkerSize',10,'MarkerFaceColor','m') ; 
hold on 
plot(xf_high_newton,0,'om','MarkerSize',10,'MarkerFaceColor','m') ; 
hold off

% check values are very high because the roots are right near poles of the
% function, indicates I have found the roots correctly
check_low = fx(xf_low_newton) 

check_high = fx(xf_high_newton)

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

function [xf_low,xf_high] = bisection(xref,eps)

    % pick bracket [a,b] s.t. a < 0 and b > 0 for xf_high 
    xah = xref - eps  ;
    xbh = xref + eps  ;
    yah = fx(xah)  ;
    ybh = fx(xbh)  ;
    
    % if the brackets are the same sign, move them until a < 0 and b > 0 
    if yah < 0 && ybh < 0 
        while ybh < 0
            xbh = xbh + eps ; 
            ybh = fx(xbh) ; 
        end
    elseif ybh > 0 && yah > 0 
        while yah > 0
            xah = xah - eps ; 
            yah = fx(xah) ; 
        end
    end
    
    % pick bracket [a,b] s.t. a > 0 and b < 0 for xf_low 
    xal = xref - eps  ;
    xbl = xref + eps  ;
    yal = fx(xal) ;
    ybl = fx(xbl) ; 
    
    % if the brackets are the same sign, move them until a > 0 and b < 0 
    if yal < 0 && ybl < 0 
        while yal < 0
            xal = xal - eps ; 
            yal = fx(xal) ; 
        end
    elseif ybl > 0 && yal > 0 
        while ybl > 0
            xbl = xbl + eps ; 
            ybl = fx(xbl) ; 
        end
    end
    
    while abs(xah - xbh) > eps
          
        % calculate midpoint (a + b)/2
        xmidh = (xah + xbh) / 2 ; 
        ymidh = fx(xmidh) ; 

        % assign midpoint to a or b depending on sign
        if ymidh > 0 && yah > 0
            xah = xmidh ; 
        elseif ymidh > 0 && ybh > 0 
            xbh = xmidh ; 
        elseif ymidh < 0 && yah < 0
            xah = xmidh ; 
        elseif ymidh < 0 && ybh < 0
            xbh = xmidh ; 
        end
        
    end
    
    xf_high = (xah + xbh) / 2 ; 
    
    while abs(xal - xbl) > eps
          
        % calculate midpoint (a + b)/2
        xmidl = (xal + xbl) / 2 ; 
        ymidl = fx(xmidl) ; 
        
        % assign midpoint to a or b depending on sign
        if ymidl > 0 && yal > 0
            xal = xmidl ; 
        elseif ymidl > 0 && ybl > 0 
            xbl = xmidl ; 
        elseif ymidl < 0 && yal < 0
            xal = xmidl ; 
        elseif ymidl < 0 && ybl < 0
            xbl = xmidl ; 
        end
        
    end
    
    xf_low = (xal + xbl) / 2 ; 
end
