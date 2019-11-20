% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3, Resubmit 1
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
iter = 5 ; 
threshDist = 2.0 ; 
inlierRatio = 0.5 ; 

[A,B,C,D] = ransac(data, minn, iter, threshDist, inlierRatio) ; 
    
% plane
[x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
z = -1*(A.*x + B.*y + D) ./ C ;  

% average distance of a point in the data set to fitted plane
d = DistPointPlane(data,A,B,C,D) ; 
E = mean(d) 

fig1 ;
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4a') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

%% part b

fid2 = fopen('cluttered_table.txt') ; 
data2 = textscan(fid2, '%10f %10f %10f') ; 
xi = data2{1,1} ; 
yi = data2{1,2} ; 
zi = data2{1,3} ; 
data2 = [xi, yi, zi] ; 

minn = 3 ; 
iter = 5 ; 
threshDist = 2.0 ; 
inlierRatio = 0.5 ; 

[A,B,C,D] = ransac(data2, minn, iter, threshDist, inlierRatio) ; 

% plane
[x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
z = -1*(A.*x + B.*y + D) ./ C ; 

figure(2)
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4b') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

fclose(fid2) ; 

%% part c
clear all ; close all ; clc ; 

fid2 = fopen('cluttered_table.txt') ; 
data2 = textscan(fid2, '%10f %10f %10f') ; 
xi = data2{1,1} ; 
yi = data2{1,2} ; 
zi = data2{1,3} ; 
data2 = [xi, yi, zi] ; 

minn = 3 ; 
iter = 10 ; 
threshDist = 1.9 ; 
inlierRatio = 0.5 ; 

[A,B,C,D] = ransac(data2, minn, iter, threshDist, inlierRatio) ; 

% plane
[x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
z = -1*(A.*x + B.*y + D) ./ C ; 

figure(2)
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4c') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

fclose(fid2) ; 

%% part d
clear all ; close all ; clc ; 

fid3 = fopen('clean_hallway.txt') ; 

data = textscan(fid3, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 
data3 = [xi, yi, zi] ; 

figure(4)

n = length(xi) / 4 ; 

for i = 1:4
    start = ((i-1)*n + 1) ; 
    finish = i*n ; 
    
    newx = xi(start:finish) ; 
    newy = yi(start:finish) ; 
    newz = zi(start:finish) ; 
    
    plot3(newx,newy,newz,'ob') ; 
    hold on 
    
    data = [newx, newy, newz] ; 
    
    minn = 3 ; 
    iter = 5 ; 
    threshDist = 3.5 ; 
    inlierRatio = 0.1 ; 

    [A,B,C,D] = ransac(data, minn, iter, threshDist, inlierRatio) ; 

    xmin = min(newx) ; 
    xmax = max(newx) ; 
    ymin = min(newy) ; 
    ymax = max(newy) ; 

    % plane
    [x, y] = meshgrid(xmin:0.01:xmax,ymin:0.001:ymax) ; 
    z = -1*(A.*x + B.*y + D) ./ C ; 

    figure(4)
    plot3(x,y,z)
    hold on
end

figure(4)
title('Problem 4d') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 


%% part e
clear all ; close all ; clc ; 

fid4 = fopen('cluttered_hallway.txt') ; 

data = textscan(fid4, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 
data4 = [xi, yi, zi] ; 

figure(5)

n = (length(xi) / 4) - 1000 ; 
E = zeros(1,4) ; 

for i = 1:4
    start = ((i-1)*n + 1) + (i*1000) ; 
    finish = i*n + (i*1000) ; 
    
    newx = xi(start:finish) ; 
    newy = yi(start:finish) ; 
    newz = zi(start:finish) ; 
    
    plot3(newx,newy,newz,'ob') ; 
    hold on 
    
    data = [newx, newy, newz] ; 
    
    minn = 3 ; 
    iter = 5 ; 
    threshDist = 3.5 ; 
    inlierRatio = 0.1 ; 

    [A,B,C,D] = ransac(data, minn, iter, threshDist, inlierRatio) ; 

    xmin = min(newx) ; 
    xmax = max(newx) ; 
    ymin = min(newy) ; 
    ymax = max(newy) ; 

    % plane
    [x, y] = meshgrid(xmin:0.01:xmax,ymin:0.001:ymax) ; 
    z = -1*(A.*x + B.*y + D) ./ C ; 

    figure(5)
    plot3(x,y,z)
    hold on
    
    % average distance of a point in the data set to fitted plane - to be
    % used as a smoothness metric
    d = DistPointPlane(data,A,B,C,D) ; 
    E(i) = mean(d)  ; 
end

% higher values for E denote rougher wall surfaces
E

figure(5)
title('Problem 4e') ; xlabel('x') ; ylabel('y') ; zlabel('z') ;

%% functions

function d = DistPointLine(xi, yi, A, B, C) 
    n = length(xi) ; 
    d = zeros(n,1) ; 
    
    for i = 1:n
        d(i) = abs(A*xi(i) + B*yi(i) + C) / sqrt(A^2 + B^2) ; 
    end
    
end

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
    
    D = var(1) ; 
    C = var(2) ; 
    B = var(3) ; 
    A = var(4) ; 
end

function [A,B,C,D] = ransac(data, minn, iter, threshDist, inlierRatio)    
    n = size(data,1) ; % total num of points
    
    bestIn = 0 ; % best fit plane with largest number of inliers
    
    % parameters for best fit plane
    A = 1 ; 
    B = 1 ; 
    C = 1 ; 
    D = 1 ; 
    
    for i = 1:iter
        
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
        end
    end
end