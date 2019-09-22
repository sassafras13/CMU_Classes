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
    
    f0 = p(1) ; 
    f1 = p(2) ; 
    f2 = p(3) ; 

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
    if b > 0 
        z = (-2*c) / (b + sqrt((b^2) - 4*a*c)) ; 
    elseif b < 0 
        z = (-2*c) / (b - sqrt((b^2) - 4*a*c)) ;
    end

    % select the root with the smallest absolute value (DO I NEED THIS?)    

    % calculate p3 
    p3 = p2 + z ; 

    % update the points for the next iteration by keeping the 2 points closest
    % to p3 and discarding the furthest point
    d = zeros(n,1) ; 
    for i = 1:length(d) 
        d(i) = abs(p3 - p(i)) ; 
    end
    
    [~,I] = min(d) ; 
    
    

    % repeat until abs(p3_new - p3_old) < epsilon
    
end
