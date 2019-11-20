% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 6
% Problem 1
% References: 
% Sommer, Pascal. "A gentle introduction to the convex hull problem." 
% https://medium.com/@pascal.sommer.ch/a-gentle-introduction-to-the-convex-hull-problem-62dfcabee90c
% Visited 11/10/2019.

%%
clc ; clear all ; close all ; 

%% 

% run a function to generate a set of n random points in 2D space
n = 10 ; 
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

% now take 3 points in a row out of this sorted stack
xconvex = zeros(length(x),2) ; 

for i = 1:length(x)

    if i == length(x) - 1
        a = x(i,1:2) ; 
        b = x(i+1,1:2) ; 
        c = x(1,1:2) ; 
    elseif i == length(x) 
        a = x(i,1:2) ; 
        b = x(1,1:2) ;
        c = x(2,1:2) ; 
    else
        % calculate the cross-product between the two vectors that make up the
        % corner 
        a = x(i,1:2) ; 
        b = x(i+1,1:2) ; 
        c = x(i+2,1:2) ; 
    end
    
    convex = isconvex(a,b,c) ; 
    
    % if the cross-product is positive then the middle point belongs on the
    % border of the convex hull and we move on to the next point
    if convex == 1
        xconvex(i,1:2) = b ; 
    elseif convex == 0
        xconvex(i,1) = NaN ; 
        xconvex(i,2) = NaN ; 
    end
end

n = length(xconvex) ; 
finalx = [] ; 

for i = 1:n
    if isnan(xconvex(i,1)) == 0
        finalx = [finalx;xconvex(i,1:2)] ; 
    end
end
 
figure(1)
plot(finalx(:,1), finalx(:,2), 'og') 
hold on
for i = 1:(length(finalx)-1)
    plot(finalx(i:i+1, 1), finalx(i:i+1, 2), '-g')
end
hold on
plot([finalx(end,1),finalx(1,1)], [finalx(end,2),finalx(1,2)],'-g')
hold on

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
    bc = [c - b, 0] ; 
    
    % calculate the cross product 
    cp = cross(ab, bc) ;
    cp = cp(3) ; 
    
    % check if cross product is > 0 and return 1 if true
    % when cross product is nonzero, vectors are NOT parallel
    if cp > 0 
        convex = 1 ; 
    elseif cp <= 0 
        convex = 0 ; 
    end
end