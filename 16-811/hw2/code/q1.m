% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 2
% Problem 1

%% Part (a) and (b)
clc ; clear all ; close all ; 

% initialize x array
x = [0, 1/8, 1/4, 1/2, 3/4, 1] ; 

% initialize fx array
fx = q1b(x) ; 

% run function to find interpolated value of fx given xdes
xdes = 1/3 ; 
pn = DivDiff(x,fx,xdes) 

% compare to calculated value 
fxcalc = q1b(xdes)

%% Part (c)
% clc ; clear all ; close all ; 

% create array of sample sizes n 
n = [2,4,40] ; 

k = length(n) ; 
figure(1)

for i = 1:k
    
    % initialize x array
    x = q1c_x(n(i)) ; 

    % initialize fx array
    fx = q1c_fx(x) ; 

    % run function to find interpolated value of fx given xdes
    xdes = 0.05 ; 
    pn = DivDiff(x,fx,xdes) 

    % compare to calculated value 
    fxcalc = q1c_fx(xdes)
    
    % plot calculated and interpolated values
    plot(n(i),pn,'or') 
    hold on
    plot(n(i),fxcalc,'xb')
    hold on
end

figure(1)
title('Interpolated vs Calculated Values (red circles = interp, blue cross = calc)') ; 
xlabel('# of Samples') ; ylabel('Value') ; 

%% Part (d)
% clc ; clear all ; close all ; 

% initialize sample size 
n = [2,4,6,8,10,12,14,16,18,20,40] ; 

% calculate maximum error for given sample size
En = MaxError(n)

figure(2)
plot(n,En,'or') 
title('Maximum Error as a Function of # of Samples') ; 
xlabel('# of Samples') ; ylabel('Maximum Error') ; 

%% Functions

% fx for Problem 1(b)
function fx = q1b(x) 
    n = length(x) ; 
    
    fx = zeros(n,1) ; 
    
    for i = 1:n
        fx(i,1) = exp(-x(i)) ; 
    end
end

% x for Problem 1(c)
function x = q1c_x(n)
    i = zeros(n,1) ; 
    x = zeros(n,1) ; 
    
    for j = 1:n
        i(j) = j-1 ; 
    end
    
    for j = 1:n
        x(j) = i(j)*(2/n) - 1 ; 
    end
end

% fx for Problem 1(c)
function fx = q1c_fx(x) 
    n = length(x) ; 
    fx = zeros(n,1) ; 
    
    for i = 1:n
        fx(i) = 1 / (1 + 16*(x(i)^2)) ; 
    end
end

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

% error function
function En = MaxError(n)

    % create the maximum error array
    k = length(n) ; 
    En = zeros(k,1) ; 
    
    % calculate maximum error for each sample size
    for i = 1:k
        
        % create x (interval size) based on i-th sample size
        x = -1:(2/n(i)):1 ; 
        
        % create error array to hold errors for every point in i-th sample
        % size
        m = length(x) ; 
        E = zeros(m,1) ; 
        
        % calculate fx array for x array
        fxcalc = q1c_fx(x) ; 
        
        % calculate error for every point in i-th sample size
        for j = 1:m
            
            % calculate pn value for current x value 
            pn_value = DivDiff(x,fxcalc,x(j)) ; 
            
            % calculate error for current x value
            E(j) = abs(fxcalc(j) - pn_value) ; 
        end
        
        En(i) = max(E) ; 
    end

end
