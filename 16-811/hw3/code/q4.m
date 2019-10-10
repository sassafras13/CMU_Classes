% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3
% Problem 4

%% 
clear all ; close all ; clc ; 

%% part a 
fid = fopen('clear_table.txt') ; 

[x, y, z, xi, yi, zi,xbar] = LSplane(fid) ;

E = LSerror(xi,yi,zi,xbar)

figure(1)
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4a') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

%% part b

fid2 = fopen('cluttered_table.txt') ; 

[x, y, z, xi, yi, zi,xbar] = LSplane(fid2) ;

E = LSerror(xi,yi,zi,xbar)

figure(2)
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4b') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

%% part c

fid2 = fopen('cluttered_table.txt') ; 

data = textscan(fid2, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 

% figure(4)
% plot(xi,yi,'ob') 
% hold on

xbar = LinearInterp(xi,yi) ; 

% x = -1.5:0.1:1.5 ; 
% fx = xbar(1) + xbar(2)*x ; 

% figure(4)
% plot(x, fx, '-r') ; 

A = xbar(2) ; % slope of the line
B = -1 ; % coefficient of y
C = xbar(1) ; % y-intercept of line

d = DistPointLine(xi, yi, A, B, C) ;

epsD = 0.025 ; % cutoff value for sorting points that are too far away

for i = 1:length(xi)
    if d(i) > epsD
        xi(i) = 0 ; 
        yi(i) = 0 ; 
        zi(i) = 0 ;
    end
end

newx = [] ; 
newy = [] ; 
newz = [] ; 

for i = 1:length(xi) 
    if xi(i) ~= 0 && yi(i) ~= 0 && zi(i) ~= 0
        newx = [newx, xi(i)] ; 
        newy = [newy, yi(i)] ; 
        newz = [newz, zi(i)] ; 
    end
end

newx = newx' ; 
newy = newy' ; 
newz = newz' ; 

xmin = -1 ; 
xmax = 1 ; 
ymin = 0.2 ; 
ymax = 0.5 ; 

[x, y, z, xbar] = LSplane_filtered(newx, newy, newz,xmin,xmax,ymin,ymax) ;

E = LSerror(newx,newy,newz,xbar)

figure(3)
% plot3(xi,yi,zi,'og') ; 
hold on
plot3(newx,newy,newz,'ob') ; 
hold on 
plot3(x,y,z) ; 
title('Problem 4c') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

%% part d

fid3 = fopen('clean_hallway.txt') ; 

data = textscan(fid3, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 

figure(4)
plot3(xi(1:4131), yi(1:4131), zi(1:4131), 'ob') ; 


n = length(xi) / 4 ; 

for i = 1:4
    start = ((i-1)*n + 1) ; 
    finish = i*n ; 
    
    newx = xi(start:finish) ; 
    newy = yi(start:finish) ; 
    newz = zi(start:finish) ; 
    
    plot3(newx,newy,newz,'ob') ; 
    hold on 
    
    xmin = min(newx) ; 
    xmax = max(newx) ; 
    ymin = min(newy) ; 
    ymax = max(newy) ; 

    [x, y, z, xbar] = LSplane_filtered(newx, newy, newz,xmin,xmax,ymin,ymax) ;
    
    figure(4)
%     plot3(x,y,z)
    f = @(x,y) xbar(1) + xbar(2)*x + xbar(3)*y ; 
    fsurf(f)
    hold on
end

figure(4)
title('Problem 4d') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 


%% part e

fid4 = fopen('cluttered_hallway.txt') ; 

data = textscan(fid4, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 

ntest = 500 ; 
xtest = xi(1:ntest) ; 
ytest = xi(1:ntest) ; 

xbar = LinearInterp(xtest,ytest) ; 

x = min(xtest):0.1:max(xtest) ; 
fx = xbar(1) + xbar(2)*x ; 

figure(5)
% plot(x, fx, '-r') ; 
% hold on
% plot(xtest,ytest,'ob') ; 
% hold on

A = xbar(2) ; % slope of the line
B = -1 ; % coefficient of y
C = xbar(1) ; % y-intercept of line

d = DistPointLine(xi, yi, A, B, C) ;

epsD = 0.985 ; % cutoff value for sorting points that are too far away

index = [0, 0, 0, 0] ; % indices of where one wall ends, another starts

j =  1 ; 
start = 6 ; 

while j < 5
    for i = start:length(xi)
        if d(i-5:i) > epsD
            index(j) = i ; 
            j = j+1 ; 
            start = i+200 ; 
            break
        end
    end
end

n = length(xi) / 4 ; 

for i = 1:4
    start = ((i-1)*n + 1) ; 
    finish = i*n ; 
    
    newx = xi(start:finish) ; 
    newy = yi(start:finish) ; 
    newz = zi(start:finish) ; 
    
    figure(5)
    plot3(newx,newy,newz,'ob') ; 
    hold on 
    
    xmin = min(newx) ; 
    xmax = max(newx) ; 
    ymin = min(newy) ; 
    ymax = max(newy) ; 

    [x, y, z, xbar] = LSplane_filtered(newx, newy, newz,xmin,xmax,ymin,ymax) ;
    
    figure(5)
    f = @(x,y) xbar(1) + xbar(2)*x + xbar(3)*y ; 
    fsurf(f)
    hold on
end

figure(5)
title('Problem 4e') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 


%% functions

function [x, y, z, xi, yi, zi,xbar] = LSplane(fid) 
    data = textscan(fid, '%10f %10f %10f') ; 
    xi = data{1,1} ; 
    yi = data{1,2} ; 
    zi = data{1,3} ; 

    phi1 = ones(length(xi),1) ; 
    phi2 = xi ; 
    phi3 = yi ; 

    A = [phi1, phi2, phi3] ; 

    % check matrix
%     invA = inv(A) 
%     nullA = null(A) 

    [U,S,V] = svd(A) ;

    % calculate inverse element by element
    invS = zeros(size(S,1),size(S,2)) ; 

    for i = 1:size(S,1)
        for j = 1:size(S,2) 
                if S(i,j) > eps
                invS(i,j) = 1/S(i,j) ; 
            else 
                invS(i,j) = 0 ; 
                end
        end
    end

    xbar = V*invS'*U*zi ;

    [x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
    
    z = xbar(1) + xbar(2)*x + xbar(3)*y ; 

end

function E = LSerror(xi,yi,zi,xbar)
    z = xbar(1) + xbar(2)*xi + xbar(3)*yi ; 
    
    p = zi - z ; 
    p2 = p.^2 ; 
    E = sum(p2) ; 
end

function [x, y, z, xbar] = LSplane_filtered(xi,yi,zi,xmin,xmax,ymin,ymax) 
    phi1 = ones(length(xi),1) ; 
    phi2 = xi ; 
    phi3 = yi ; 

    A = [phi1, phi2, phi3] ; 

    % check matrix
%     invA = inv(A) 
%     nullA = null(A) 

    [U,S,V] = svd(A) ;

    % calculate inverse element by element
    invS = zeros(size(S,1),size(S,2)) ; 

    for i = 1:size(S,1)
        for j = 1:size(S,2) 
                if S(i,j) > eps
                invS(i,j) = 1/S(i,j) ; 
            else 
                invS(i,j) = 0 ; 
                end
        end
    end

    xbar = V*invS'*U*zi ;

    [x, y] = meshgrid(xmin:0.01:xmax,ymin:0.001:ymax) ; 
    
    z = xbar(1) + xbar(2)*x + xbar(3)*y ; 

end

function xbar = LinearInterp(xi,yi)

    phi1 = ones(length(yi), 1) ; 
    phi2 = xi ; 

    A = [phi1, phi2] ; 

    % check matrix
    % invA = inv(A) 
    % nullA = null(A) 

    [U,S,V] = svd(A) ;

    % calculate inverse element by element
    invS = zeros(size(S,1),size(S,2)) ; 

    for i = 1:size(S,1)
        for j = 1:size(S,2) 
                if S(i,j) > eps
                invS(i,j) = 1/S(i,j) ; 
            else 
                invS(i,j) = 0 ; 
                end
        end
    end

    xbar = V*invS'*U*yi ; 

end

function d = DistPointLine(xi, yi, A, B, C) 
    n = length(xi) ; 
    d = zeros(n,1) ; 
    
    for i = 1:n
        d(i) = abs(A*xi(i) + B*yi(i) + C) / sqrt(A^2 + B^2) ; 
    end
    
end
