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
x0 = 1 ; 
y0 = 1 ; 
p0 = [x0,y0] ; 

% pick 3 paths such that the given starting point lies inside a triangle formed
% by the starting points of these 3 paths

% build a new path p by performing a weighted sum of the 3 paths. The
% weights should be constant with time. P(0) should just be the given
% starting point. 

% indices of the chosen paths (+1 for y-coordinates)
i = 1 ; 
j = 3 ; 
k = 5 ; 

% weights for each of the chosen paths
a1 = 0.5 ; 
a2 = 0.25 ; 
a3 = 0.25 ; 

% calculated path
px = [x0, a1*data(i,2:end) + a2*data(j,2:end) + a3*data(k,2:end)] ; 
py = [y0, a1*data(i+1,2:end) + a2*data(j+1,2:end) + a3*data(k+1,2:end)] ; 

% use interpolation to find the values for p(t) at some given time scale
% (we will need to justify why we picked the time scale and the
% interpolation method)
tbaseline = 1:1:50 ; 
dt = 0.1 ; 
t = 1:dt:50 ; 

% I will use MATLAB's interp1 (linear interpolation)
px_interp = interp1(tbaseline,px,t) ; 
py_interp = interp1(tbaseline,py,t) ; 

% p(t) should never enter the ring of fire

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
legend('Ring of Fire','Path 1','Path 2','Path 3','Interpolated Path','Location','southeast') ; 
xlim([0 12]) ; ylim([0 12]) ; 

%% Functions

% divided difference function
function pn = DivDiff(x,fx,xdes)

    % get number of sample points
    n = length(x) ; 
    ncheck = length(fx) ; 
    
    if n ~= ncheck
        f = msgbox("Sample point arrays x and fx do not have the same length.","Error") ; 
        return ; 
    end
    
    % calculate b
    % the structure of b is: 
    % b = [fx ,  0 ,   0 ,   0 ...
    %      fx2 , b11 , 0 ,   0 ...
    %      fx3 , b12 , b21 , 0 ...
    %      fx4 , b13 , b22 , b31 ...]
    % b = j rows x i columns
    
    % initialize b
    b = zeros(n,n) ; 
    
    % first column of b uses fx values
    for j = 1:n
            b(j,1) = fx(j) ; 
    end
    
    % columns 2:n are calculated values of b
    for i = 2:n
        for j = (1 + (i-1)):n 
            b(j,i) = ( b(j,(i-1)) - b((j-1),(i-1)) ) / ( x(j) - x(j-1) ) ; 
        end
    end        
    
    % find starting x values
    % search for xlow which is less than xdes and xhi which is greater than
    % xdes
    % xlow and xhi should bracket xdes
    % save the indices of xlow and xhi, these should be used to build the
    % polynomial starting with index of xlow, and expanding up to larger
    % values
    for j = 1:(n-1)
        xlow = x(j) ; 
        jlow = j ; 
        xhigh = x(j+1) ; 
        
        if xlow < xdes && xhigh > xdes
            break
        end
    end
    
    % build polynomial px
    % solve px for xdes
    px = b(jlow,1) ; 
    xmultiplier = 1 ; 
    
    j = jlow ; % start looping through jlow-th row and go down
    i = 2 ; 
    
    while j <= n && i <= n 
        xmultiplier = xmultiplier * (xdes - x(j)) ; % calculate (x - xi)
        px = px + b(j,i)*xmultiplier ;   % add new term to polynomial
        j = j + 1 ; % increase row index by 1 
        i = i + 1 ; % increase column index by 1

    end
    
    pn = px ; 

end

% plot a circle
function [xc,yc] = circle(x,y,r)
    theta = 0:pi/50:2*pi ; 
    xc = r*cos(theta) + x ; 
    yc = r*sin(theta) + y ; 
end