% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 5
% References: Mathews, J. and Fink, K. "Numerical Methods Using Matlab, 4th
% Edition." Prentice-Hall Inc., 2004.
% <http://mathfaculty.fullerton.edu/mathews/n2003/mullersmethod/MullersMethodProof.pdf>
% Visited 09/21/2019. 

% Velix, Oscar. "Muller's Method."
% <https://www.youtube.com/watch?v=XIIEjwtkONc> Visited 10/14/2019.

%% 
clc ; clear all ; close all ; 

%% Plot function to understand how it looks

x = -50:0.1:50 ; 
y = fx(x) ;  
line = zeros(length(x),1) ; 

f = figure(1) ; 
plot(x,y,'-b') ; 
hold on
plot(x,line,'-r'); 
hold on 
xlabel('x') ; ylabel('y') ; 

% from inspection it looks like the root is at x = 3

%% Perform Muller's method 

% create initial 3 points
p0 = -40 ; 
p1 = -30 ; 
p2 = -20 ; 

p = [p0, p1, p2] ; 

% obtain fx values for initial 3 points
f0 = fx(p0) ; 
f1 = fx(p1) ; 
f2 = fx(p2) ; 

fp = [f0, f1, f2] ; 

% call the function 
xf = mullerMethod2(p,fp)

% first root is 3, so collapse and find new roots
q0 = -10 ; 
q1 = 3 ; 
q2 = 10 ;

q = [q0, q1, q2] ; 

fq0 = qx(q0) ; 
fq1 = qx(q1) ; 
fq2 = qx(q2) ; 

fq = [fq0, fq1, fq2] ; 

xf2 = mullerMethod3(q,fq)  

%% Functions 

function y = fx(x) 
    y = x.^3 - 5.*(x.^2) + 11.*x - 15 ; 
end

function y = qx(x)
    y = x.^2 - 2.*x + 5 ; 
end

function xf = mullerMethod2(p,fp) 

    % get length of p 
    n = length(p) ; 
    
    % extract values
    p0 = p(1) ; 
    p1 = p(2) ; 
    p2 = p(3) ; 
    
    f0 = fp(1) ; 
    f1 = fp(2) ; 
    f2 = fp(3) ; 
    
    % epsilon
    eps = 1e-1 ; 
    
    % define p3new and p3old
    p3old = p2 ; 
    p3new = 0 ; 
    
    % define check
    check = abs(p3new - p3old) ; 

    % repeat until abs(p3_new - p3_old) < epsilon
    while check > eps
        
        % calculate q
        q = (p2 - p1) / (p1 - p0) ; 
        
        a = q*f2 - q*(1+q)*f1 + (q^2)*f0 ; 
        b = (2*q + 1)*f2 - ((1+q)^2)*f1 + (q^2)*f0 ; 
        c = (1 + q)*f2 ; 
        
        denom1 = b + sqrt((b^2) - 4*a*c) ; 
        denom2 = b - sqrt((b^2) - 4*a*c) ; 
        
        denom = max(denom1, denom2) ; 
        
        p3 = p2 - (p2 - p1)*((2*c)/denom) ; 
        
        % define p3new and p3old
        p3old = p2 ; 
        p3new = p3 ; 

        % define check
        check = abs(p3new - p3old) ; 
    
        % update the points for the next iteration by keeping the 2 points closest
        % to p3 and discarding the furthest point 
        pfull = [p0, p1, p2, p3] ; 
        
        d = zeros(length(pfull),1) ; 
        
        for i = 1:length(pfull) 
            d(i) = abs(p3 - pfull(i)) ; 
        end
        
        % delete furthest point from array
        [~,I] = max(d) ; 
        pfull(I) = [] ; 

        % assign the remaining points for the next iteration
        %pfull = sort(pfull) ;         
        p0 = pfull(1) ; 
        p1 = pfull(2) ; 
        p2 = pfull(3) ;
        
        f0 = fx(p0) ; 
        f1 = fx(p1) ; 
        f2 = fx(p2) ; 
    end
    
    xf = p2 ; 
            
end

function xf = mullerMethod3(p,fp) 

    % get length of p 
    n = length(p) ; 
    
    % extract values
    p0 = p(1) ; 
    p1 = p(2) ; 
    p2 = p(3) ; 
    
    f0 = fp(1) ; 
    f1 = fp(2) ; 
    f2 = fp(3) ; 
    
    % epsilon
    eps = 1e-1 ; 
    
    % define p3new and p3old
    p3old = p2 ; 
    p3new = 0 ; 
    
    % define check
    check = abs(p3new - p3old) ; 

    % repeat until abs(p3_new - p3_old) < epsilon
    while check > eps
        
        % calculate q
        q = (p2 - p1) / (p1 - p0) ; 
        
        a = q*f2 - q*(1+q)*f1 + (q^2)*f0 ; 
        b = (2*q + 1)*f2 - ((1+q)^2)*f1 + (q^2)*f0 ; 
        c = (1 + q)*f2 ; 
        
        denom1 = b + sqrt((b^2) - 4*a*c) ; 
        denom2 = b - sqrt((b^2) - 4*a*c) ; 
        
        denom = max(denom1, denom2) ; 
        
        p3 = p2 - (p2 - p1)*((2*c)/denom) ; 
        
        % define p3new and p3old
        p3old = p2 ; 
        p3new = p3 ; 

        % define check
        check = abs(p3new - p3old) ; 
    
        % update the points for the next iteration by keeping the 2 points closest
        % to p3 and discarding the furthest point 
        pfull = [p0, p1, p2, p3] ; 
        
        d = zeros(length(pfull),1) ; 
        
        for i = 1:length(pfull) 
            d(i) = abs(p3 - pfull(i)) ; 
        end
        
        % delete furthest point from array
        [~,I] = max(d) ; 
        pfull(I) = [] ; 

        % assign the remaining points for the next iteration
        %pfull = sort(pfull) ;         
        p0 = pfull(1) ; 
        p1 = pfull(2) ; 
        p2 = pfull(3) ;
        
        f0 = qx(p0) ; 
        f1 = qx(p1) ; 
        f2 = qx(p2) ; 
    end
    
    xf = p2 ; 
            
end
