% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3, Resubmit 2
% Problem 4
% References: 
% (1) https://en.wikipedia.org/wiki/Random_sample_consensus

%% 
clear all ; close all ; clc ; 

%% part a 
fig1 = figure(1) ;

fid = fopen('clear_table.txt') ; 
data = textscan(fid, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 
data = [xi, yi, zi] ; 

minn = 3 ; 
iter = 100 ; 
threshDist = 2.0 ; 
inlierRatio = 0.5 ; 

[A,B,C,D] = ransac(data, minn, iter, threshDist, inlierRatio) ; 
    
% plane
[x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
z = -1*(A.*x + B.*y + D) ./ C ;  

% average distance of a point in the data set to fitted plane
d = DistPointPlane(data,A,B,C,D) ; 
E1 = mean(d) 

fig1 ;
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4a') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

%% part e
clear all ; close all ; clc ; 

global data

eps = 0.1 ; % cutoff for one plane

fid4 = fopen('cluttered_hallway.txt') ; 

data = textscan(fid4, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 
data = [xi, yi, zi] ; 

xmin = -1.5 ; 
xmax = 1.5 ; 
ymin = -0.5 ; 
ymax = 1.5 ; 
zmin = 1 ; 
zmax = 6 ; 

% figure(5)
% plot3(xi,yi,zi,'ob') ; 
% hold on 

minn = 3 ; 
iter = 200 ; 
threshDist = 1.1 ; 
inlierRatio = 0.5 ; 

for i = 1:4
    [A,B,C,D] = ransac(minn, iter, threshDist, inlierRatio) ;

    [x, y] = meshgrid(xmin:0.01:xmax,ymin:0.001:ymax) ; 
    z = -1*(A.*x + B.*y + D) ./ C ; 

    figure(5)
    plot3(x,y,z)
    hold on
end

figure(5)
title('Problem 4e') ; xlabel('x') ; ylabel('y') ; zlabel('z') ;

%% functions

function d = DistPointPlane(sample,A,B,C,D)
    d = abs(A*sample(:,1) + B*sample(:,2) + C*sample(:,3) + D) / sqrt(A^2 + B^2 + C^2) ; 
end

function [A,B,C,D] = fitPlane(sample) 
    P1 = sample(1,:) ; 
    P2 = sample(2,:) ; 
    P3 = sample(3,:) ; 
    
    normal = cross(P1-P2, P1-P3) ; 
    
    syms x y z 'real'
    P = [x,y,z] ; 
    realdot = @(u,v) u*transpose(v) ; 
    planefunction = realdot(normal, P-P1) ; 
    
    var = double(coeffs(planefunction)) ; 
    
    if size(var,2) == 1
        var = [var 0 0 0] ; 
    elseif size(var,2) == 2 
        var = [var 0 0] ; 
    elseif size(var,2) == 3
        var = [var 0] ; 
    elseif size(var,2) == 0
        var = [0 0 0 0] ;
    end
    
    D = var(1) ; 
    C = var(2) ; 
    B = var(3) ; 
    A = var(4) ; 
end

function [A,B,C,D] = ransac(minn, iter, threshDist, inlierRatio)        
    global data
    bestIn = 0 ; % best fit plane with largest number of inliers
    
    % parameters for best fit plane
    A = 1 ; 
    B = 1 ; 
    C = 1 ; 
    D = 1 ; 
    
    for i = 1:iter
        n = size(data,1) ; % total num of points
        
        % randomly pick 3 points
        index = randperm(n,minn) ; 
        sample = data(index,:) ; 
        
        % calculate the distances between all points and the fit plane
        d = DistPointPlane(data,A,B,C,D) ; 
        
        % plane
        [x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
        z = -1*(A.*x + B.*y + D) ./ C ; 
        
        % compute inliers with distances smaller than the threshold
        inlierIndex = find(abs(d) <= threshDist) ; 
        inlierNum = length(inlierIndex) ; 
        
        % update the number of inliers and fitting model if a better model
        % is found
        if inlierNum >= round(inlierRatio*n) && inlierNum >= bestIn
            bestIn = inlierNum ; 
            [A,B,C,D] = fitPlane(sample) ;
             
            
%             figure(5) 
%             plot3(data(inlierIndex,1), data(inlierIndex,2), data(inlierIndex,3),'ob')
%             hold on
            
            data(index,:) = [] ;
        end
    end
end