% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 6
% Problem 1
% References: 
% Sommer, Pascal. "A gentle introduction to the convex hull problem." 
% https://medium.com/@pascal.sommer.ch/a-gentle-introduction-to-the-convex-hull-problem-62dfcabee90c
% Visited 11/10/2019.

% Wikipedia. "Graham scan." https://en.wikipedia.org/wiki/Graham_scan
% Visited 11/20/2019.

%%
clc ; clear all ; close all ; 

%% 

% run a function to generate a set of n random points in 2D space
n = 100000 ; 
x = randpoints(n) ; 

% plot the points
figure(1)
plot(x(:,1), x(:,2), 'ob')
hold on
xlabel('x') ; ylabel('y') ; 

% look for the point that has the smallest y value
% if there are multiple such points, choose the one with the largest x
% value
[originy, ii] = min(x(:,2)) ; % ii is the index of the minimum value

% plot the minimum value in red
figure(1)
plot(x(ii,1), x(ii,2), 'or') 
hold on

% set the minimum value as the origin and remove it from the list of points
origin = [x(ii,1), originy] ; 
x(ii,:) = [] ; 

% calculate the angle between the origin and all other points in the
% set
angles = calcangles(origin, x) ; 

% sort all the points by the size of the angles calculated above
% if 2 points have the same angle we choose the one closer to the chosen
% point
x = [x, angles] ; % adds angles data as a column to x matrix of points
x = sortrows(x, 3) ; % sorts matrix x by angle size
x = [origin, 0 ; x] ; % add back in a row for the origin at the beginning 

% x =     [7.7102    8.9504         0 ;
%    37.5433   12.9271    0.1325;
%    83.1433   43.7616    0.4324;
%    62.0553   41.3396    0.5375;
%    97.8770   65.8992    0.5633;
%    70.2556   56.5651    0.6507;
%    46.2759   42.2690    0.7125;
%    95.2576   96.9590    0.7880;
%    92.0795   94.8999    0.7947;
%    12.0622   64.6387    1.4928 ] ; 
% 
% figure(1)
% plot(x(:,1),x(:,2),'ob') ; 
% hold on 

% now take 3 points in a row out of this sorted stack
stack = [] ;
stack = [stack ; x(1:2,:)] ;

for i = 3:length(x) 
    while length(stack) > 1 && isconvex(stack(end-1,1:2),stack(end,1:2),x(i,1:2)) == 0 
        stack(end,:) = [] ; 
    end
    stack = [stack ; x(i,:)] ; 
end

finalx = stack ; 

figure(1)
plot(finalx(1:end-1,1), finalx(1:end-1,2), 'or') 
hold on
for i = 1:(length(finalx)-1)
    plot(finalx(i:i+1, 1), finalx(i:i+1, 2), '-r','LineWidth',3)
end
hold on
plot([finalx(end,1),finalx(1,1)], [finalx(end,2),finalx(1,2)],'-r','LineWidth',3)
hold on

figure(1)
title('HW6 Q1') ; xlabel('X') ; ylabel('Y') ; 

%% Functions

function x = randpoints(n)
    x = 100*rand([n,2]) ; % scale the random numbers so they are between 0 and 100
end

function angles = calcangles(origin, x) 
    angles = zeros(length(x),1) ; 
    horz = [1,0] ; 
    
    for i = 1:length(x)
        
        % create a vector using the vector to the point x
        vec = x(i,:) - origin ; 
        
        % calculate the angle between the two vectors in radians
        angles(i) = acos( dot(horz,vec) / (norm(horz) * norm(vec)) ) ; 
    end
end

function convex = isconvex(a,b,c)
    % build the vectors
    ab = [b - a, 0] ; 
    ac = [c - a, 0] ; 
    
    % calculate the cross product 
    cp = cross(ab, ac) ;
    cp = cp(3) ; 
    
    % check if cross product is > 0 and return 1 if true
    % when cross product is nonzero, vectors are NOT parallel
    if cp > 0 
        convex = 1 ; 
    elseif cp <= 0 
        convex = 0 ; 
    end
end