% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 8
% References: 
% [1]
% https://www.mathworks.com/matlabcentral/answers/98665-how-do-i-plot-a-circle-with-a-given-radius-and-center
% Visited 09/22/2019

%% 
clc ; clear all ; close all ; 

%% Setup 

% import all the paths from the provided .txt file
data = importdata('paths.txt') ; 

% the structure of data is: 
% data = [x1 ... (i = 1) 
%         y1 ... (i = 2) 
%         x2 ... (i = 3) 
%         y2 ... (i = 4)

% if I want the i-th path then I want the (2*i - 1) row for x values and
% the (2*i) row for y values

%% Algorithm

% take a starting location for t = 0 
% p0 = [0.8, 1.8] ; 
% p0 = [2.2, 1.0] ; 
p0 = [2.7, 1.4] ; 

% pick 3 paths such that the given starting point lies inside a triangle formed
% by the starting points of these 3 paths

% the 3 paths be within a distance eps of each other
eps = 3 ; 

% indices of the chosen paths (+1 for y-coordinates)
i = 1 ; 
j = 3 ; 
k = 5 ; 

% build x array and y array
xarray = [data(i,1), data(j,1), data(k,1)] ; 
yarray = [data(i+1,1), data(j+1,1), data(k+1,1)] ; 

% check that the starting points of the 3 paths are within eps distance of each other
xcheck = max(xarray) - min(xarray) ; 
ycheck = max(yarray) - min(yarray) ; 

% check that all 3 paths are either above or below the ring of fire
y1 = interp1(data(i,:),data(i+1,:),5) ; 
y2 = interp1(data(j,:),data(j+1,:),5) ; 
y3 = interp1(data(k,:),data(k+1,:),5) ; 

% p(t) should never enter the ring of fire
if (y1 > 7)  && (y2 > 7) && (y3 > 7)
    circleCheck = 1 ; 
elseif (y1 < 3) && (y2 < 3) && (y3 < 3)
    circleCheck = 1 ; 
else 
    circleCheck = 0 ; 
end

result = triangleCheck(xarray, yarray, p0) ;  

m = 1 ; 

% sweep through all points 
it = 1:2:97 ; 
jt = 1:2:97 ; 
kt = 1:2:97 ; 

[I,J,K] = meshgrid(it,jt,kt) ; 
I = reshape(I,[117649,1]) ; 
J = reshape(J,[117649,1]) ; 
K = reshape(K,[117649,1]) ; 

while (result == 0) || (xcheck > eps) || (ycheck > eps) || (circleCheck == 0)
        
    i = I(m) ; 
    j = J(m) ; 
    k = K(m) ; 
    
    % build x array and y array
    xarray = [data(i,1), data(j,1), data(k,1)] ; 
    yarray = [data(i+1,1), data(j+1,1), data(k+1,1)] ; 

    % check that the points are within eps distance of each other
    xcheck = max(xarray) - min(xarray) ; 
    ycheck = max(yarray) - min(yarray) ; 
    
    % check that all 3 paths are either above or below the ring of fire
    y1 = interp1(data(i,:),data(i+1,:),5) ; 
    y2 = interp1(data(j,:),data(j+1,:),5) ; 
    y3 = interp1(data(k,:),data(k+1,:),5) ; 

    % p(t) should never enter the ring of fire
    if (y1 > 7)  && (y2 > 7) && (y3 > 7)
        circleCheck = 1 ; 
    elseif (y1 < 3) && (y2 < 3) && (y3 < 3)
        circleCheck = 1 ; 
    else 
        circleCheck = 0 ; 
    end

    result = triangleCheck(xarray, yarray, p0) ; 
    
    m = m + 1 ; 
end

% build a new path p by performing a weighted sum of the 3 paths. The
% weights should be constant with time. P(0) should just be the given
% starting point. 

% weights for each of the chosen paths
a1 = 0.5 ; 
a2 = 0.25 ; 
a3 = 0.25 ; 

% calculated path
px = [p0(1), a1*data(i,2:end) + a2*data(j,2:end) + a3*data(k,2:end)] ; 
py = [p0(2), a1*data(i+1,2:end) + a2*data(j+1,2:end) + a3*data(k+1,2:end)] ; 

% use interpolation to find the values for p(t) at some given time scale
% (we will need to justify why we picked the time scale and the
% interpolation method)
tbaseline = 1:1:50 ; 
dt = 0.1 ; 
t = 1:dt:50 ; 

% I will use MATLAB's interp1 (linear interpolation)
px_interp = interp1(tbaseline,px,t) ; 
py_interp = interp1(tbaseline,py,t) ; 


%% Demonstration

% interpolate paths for the provided starting points

% plot the ring of fire, the 3 paths that are being interpolated, as well
% as the interpolated path. show you do not intersect the ring of fire. 

% ring of fire
[xc,yc] = circle(5,5,1.5) ; 

% 3 chosen paths
px1 = data(i,:) ; 
py1 = data(i+1,:) ; 
px2 = data(j,:) ; 
py2 = data(j+1,:) ; 
px3 = data(k,:) ; 
py3 = data(k+1,:) ; 

% plot
figure(1)
plot(xc,yc,'-r') ; 
hold on 
plot(px1,py1,'--b') ; 
hold on 
plot(px2,py2,'-xb') ; 
hold on
plot(px3,py3,'-ob') ; 
hold on
plot(px_interp,py_interp,'-g') ; 
hold on 
xlabel('X') ; ylabel('Y') ; 
legend('Ring of Fire','Path 1','Path 2','Path 3','Interpolated Path','Location','northwest') ; 
xlim([0 12]) ; ylim([0 12]) ; 
title('Starting Point (2.7, 1.4)') ; 
axis equal

%% Functions

% plot a circle
function [xc,yc] = circle(x,y,r)
    theta = 0:pi/50:2*pi ; 
    xc = r*cos(theta) + x ; 
    yc = r*sin(theta) + y ; 
end

% check if a point is inside a triangle or outside
function result = triangleCheck(xarray, yarray, p0)
    A = [xarray ; 
         yarray ; 
         ones(1,3) ] ; 
     
    b = [p0, 1] ; 
    
    v = inv(A)*b' ; 
    
    if (sum(v>0) == 3) && (sum(v) == 1)
        result = 1 ; 
    else 
        result = 0 ; 
    end
end