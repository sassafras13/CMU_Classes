% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 3, Resubmit 2
% Problem 4
% References: 
% (1) https://en.wikipedia.org/wiki/Random_sample_consensus
% (2) http://www.cse.yorku.ca/~kosta/CompVis_Notes/ransac.pdf
% (3) http://www.cse.psu.edu/~rtc12/CSE486/lecture15.pdf

%% part a 
clear all ; close all ; clc ; 

fid = fopen('clear_table.txt') ; 
data = textscan(fid, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 
data = [xi, yi, zi] ; 

[x, y, z, xi, yi, zi,A,B,C,D] = LSplane(data) ;

d = DistPointPlane(data,A,B,C,D) ; 

format long
e = mean(d) 

figure(1)
plot3(xi,yi,zi,'ob') ; 
hold on
plot3(x,y,z) ; 
title('Problem 4a') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 


%% part d
clear all ; close all ; clc ; 
      
% fid2 = fopen('cluttered_table.txt') ; 
% data2 = textscan(fid2, '%10f %10f %10f') ; 
% xi = data2{1,1} ; 
% yi = data2{1,2} ; 
% zi = data2{1,3} ; 
% data2 = [xi, yi, zi] ; 

% fid3 = fopen('clean_hallway.txt') ; 
% 
% data = textscan(fid3, '%10f %10f %10f') ; 
% xi = data{1,1} ; 
% yi = data{1,2} ; 
% zi = data{1,3} ; 
% data3 = [xi, yi, zi] ; 

fid4 = fopen('cluttered_hallway.txt') ; 

data = textscan(fid4, '%10f %10f %10f') ; 
xi = data{1,1} ; 
yi = data{1,2} ; 
zi = data{1,3} ; 
data4 = [xi, yi, zi] ; 

figure(4)
plot3(xi,yi,zi,'ob')
hold on
xlabel('x') ; ylabel('y') ; 
    
minn = 3 ; 
iter = 5 ; 
% threshDist = 0.015 ; % tune
% inlierRatio = 0.21 ; % tune

Num = size(data4,1) ; 

i = 1 ; 

while i < 5

    if i < 3 
        threshDist = 0.06 ; 
        inlierRatio = 0.26 ; 
    else
        threshDist = 0.015 ; 
        inlierRatio = 0.04 ; 
    end
    
    [A,B,C,D,inlierIndex,errors,inliers] = ransac(data4, minn, iter, threshDist, inlierRatio, Num) ; 

    % plane
    xmin = min(inliers(:,1)) ; 
    xmax = max(inliers(:,1)) ; 
    ymin = min(inliers(:,2)) ;
    ymax = max(inliers(:,2)) ; 
    zmin = min(inliers(:,3)) ;
    zmax = max(inliers(:,3)) ;
    [x, y] = meshgrid(xmin:0.01:xmax,ymin:0.001:ymax) ; 
    z = -1*(A.*x + B.*y + D) ./ C ; 

    figure(4)
    plot3(x,y,z) ; 
    hold on 
    
    numInliers = size(inlierIndex,1) ;
    smoothness = errors 
    
    data4(inlierIndex,:) = [] ; 
    i = i + 1 ; 
end

figure(4)
title('Problem 4e') ; xlabel('x') ; ylabel('y') ; zlabel('z') ; 

%% functions

function [x, y, z, xi, yi, zi,A,B,C,D] = LSplane(data) 
    xi = data(:,1) ; 
    yi = data(:,2) ; 
    zi = data(:,3) ; 

    phi1 = ones(length(xi),1) ; 
    phi2 = xi ; 
    phi3 = yi ; 

    A = [phi1, phi2, phi3] ; 


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
    
    xbar = V*invS'*U'*zi ;

    [x, y] = meshgrid(-1.5:0.01:1.5,0.2:0.001:0.5) ; 
    
    z = xbar(1) + xbar(2)*x + xbar(3)*y ; 
    
    D = xbar(1) ; 
    A = xbar(2) ; 
    B = xbar(3) ; 
    C = -1 ; 
end

function d = DistPointPlane(sample,A,B,C,D)
    d = abs(A.*sample(:,1) + B.*sample(:,2) + C.*sample(:,3) + D) / sqrt(A.^2 + B.^2 + C.^2) ; 
end

function [A,B,C,D,inlierIndex,errors,inliers] = ransac(data, minn, iter, threshDist, inlierRatio, Num)        

    % save all plane coefficients
    planeCoeff = zeros(iter,4) ; 
    
    % save all errors
    errors = zeros(iter,1) ; 
    
    % save all inliers
    allInliers = cell(iter,1) ; 
    
    % save all inlier indices
    allInlierInds = cell(iter,1) ; 
    
    for i = 1:iter
        
        % create good parameter to check if # inliers > threshold
        good = 0 ; 
        j = 0 ; 
        while good == 0 && j < 1000 
            n = size(data,1) ; % total num of points
            
            % randomly pick 3 points
            index = randperm(n,minn) ; 
            sample = data(index,:) ; 

            % fit a plane to these points
            [~, ~, ~, ~, ~, ~,A,B,C,D] = LSplane(sample) ; 

            % calculate the distances between all points and the fit plane
            d = DistPointPlane(data,A,B,C,D) ; 

            % compute inliers with distances smaller than the threshold
            inlierIndex = find(abs(d) <= threshDist) ; 
            inlierNum = length(inlierIndex) ; 

            % create an array containing all the inliers
            inliers = data(inlierIndex,:) ; 
            
            j = j + 1 ;
            
            % update the number of inliers and fitting model if a better model
            % is found
            if inlierNum >= round(inlierRatio*Num) || j == 1000
                % fit the plane to all the inliers
                [~, ~, ~, ~, ~, ~,A,B,C,D] = LSplane(inliers) ;
                
                % update the check parameter as true
                good = 1  ;
                
                % calculate the fitting error 
                d = DistPointPlane(data,A,B,C,D) ; 
                e = mean(d) ;
                
                % save plane coefficients
                planeCoeff(i,:) = [A,B,C,D] ; 
                
                % save all errors
                errors(i,1) = e ; 
                
                % save all inliers
                allInliers{i} = inliers ; 
                
                % save all inlier indices
                allInlierInds{i} = inlierIndex ; 
                
            end
        end
    end
    
    [~,I] = min(errors) ; 
    coeffs = planeCoeff(I,:) ; 
    A = coeffs(1) ; 
    B = coeffs(2) ; 
    C = coeffs(3) ; 
    D = coeffs(4) ; 
    inlierIndex = allInlierInds{I} ;     
    inliers = allInliers{I} ; 
    errors = errors(I,1) ; 
end