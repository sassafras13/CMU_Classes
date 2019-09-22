% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 5
% References: Mathews, J. and Fink, K. "Numerical Methods Using Matlab, 4th
% Edition." Prentice-Hall Inc., 2004.
% <http://mathfaculty.fullerton.edu/mathews/n2003/mullersmethod/MullersMethodProof.pdf>
% Visited 09/21/2019. 

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
xf = mullerMethod(p,fp)

%% Functions 

function y = fx(x) 
    y = x.^3 - 5.*(x.^2) + 11.*x - 15 ; 
end

function xf = mullerMethod(p,fp) 

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
    eps = 1e-3 ; 
    
    % define p3new and p3old
    p3old = p2 ; 
    p3new = 0 ; 
    
    % define check
    check = abs(p3new - p3old) ; 

    % repeat until abs(p3_new - p3_old) < epsilon
    while check > eps
        
        % define h0 and h1
        h0 = p0 - p2 ; 
        h1 = p1 - p2 ; 

        % set c = f2
        c = f2 ; 

        % define e0 and e1
        e0 = f0 - c ; 
        e1 = f1 - c ; 

        % calculate a and b 
        a = (e0*h1 - e1*h0) / (h1*(h0^2) - h0*(h1^2)) ; 
        b = (e1*(h0^2) - e0*(h1^2)) / (h1*(h0^2) - h0*(h1^2)) ; 

        % find the roots using modified quadratic formula
        % if b > 0 use + ; if b < 0 use -
%         if b > 0 
%             z = (-2*c) / (b + sqrt((b^2) - 4*a*c)) ; 
%         elseif b < 0 
%             z = (-2*c) / (b - sqrt((b^2) - 4*a*c)) ;
%         end
        z1 = (-2*c) / (b + sqrt((b^2) - 4*a*c)) ;
        z2 = (-2*c) / (b - sqrt((b^2) - 4*a*c)) ;

        % select the root with the smallest absolute value 
        if imag(z1) > imag(z2)
            z = z2 ; 
        elseif imag(z2) > imag(z1)
            z = z1 ; 
        end
        
        % use the absolute value (i.e. magnitude) of the root
        z = z ; 

        % calculate p3 
        p3 = p2 + z ; 
        
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
    end
    
    xf = p2 ; 
            
end
