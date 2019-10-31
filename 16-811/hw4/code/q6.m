% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 4
% Problem 6
% NOTE: code is based on sample code provided by class instructor
% References: 
% [1] https://en.wikipedia.org/wiki/Gradient_descent

%% 
clc ; clear all ; close all ; 

%% create obstacle field
close all
clear all
waypoints=300;
N=101;
OBST = [20,30;60,40;70,85]; % defining the obstacles
epsilon = [25; 20; 30];

% creating the obstacle cost for the whole map
obs_cost = double(zeros(N));
for i=1:size(OBST,1)

    t = zeros(N);
    t(OBST(i,1),OBST(i,2)) = 1; %point obstacles
    
    t_cost = double(bwdist(t));
    t_cost(t_cost>epsilon(i))=epsilon(i);
    t_cost = 1/(2*epsilon(i))*(t_cost-epsilon(i)).^2;
    
    obs_cost = obs_cost + t_cost(1:N, 1:N);
end

% plot the obstacle map
figure(1)
imagesc(obs_cost')
% obstacle cost gradient
gx = diff(double(obs_cost),1,1);
gy = diff(double(obs_cost),1,2);
hold on


figure(1);
surface(1:N,1:N,double(obs_cost'));
xlabel('X')
hold on;

%% initial path
%world params
SX = 10; % START
SY = 10;
GX = 90; % GOAL
GY = 90;

% plot the initial trajectory as a straight line through the map
traj = zeros(2,waypoints);
traj(1,1)=SX;
traj(2,1)=SY;
dist_x = GX-SX;
dist_y = GY-SY;
for i=2:waypoints
    traj(1,i)=traj(1,i-1)+dist_x/(waypoints-1);
    traj(2,i)=traj(2,i-1)+dist_y/(waypoints-1);
end

% turn the trajectory into a 3D path by coming up with the z-height at each
% x,y position
path_init = traj';
tt=size(path_init,1);
path_init_values = zeros(size(path_init,1),1);
for i=1:tt
    path_init_values(i)=obs_cost(floor(path_init(i,1)),floor(path_init(i,2)));
end

% plot the x,y position and then the z-height
plot3(path_init(:,1),path_init(:,2),path_init_values,'.r','MarkerSize',20);
hold on

path = path_init;

%% part (a) 

gamma = 5 ; % provided in question 
maxiters = 5 ; % max number of iterations 
eps = 0.1 ; % termination tolerance

path_values = path_init_values' ; 
mse = 0.5 * sum(path_values.^2) ; 
iter = 1 ; 

% run a while loop until convergence
while iter < maxiters && mse > eps 

    % new trajectory
    newtraj = zeros(2,waypoints);
    newtraj(1,1)=SX;
    newtraj(2,1)=SY;

    for i = 2:waypoints
       % find the gradient for the current point in the trajectory
       r = floor(traj(1,i-1)) ; % x position
       c = floor(traj(2,i-1)) ; % y position

       % I am optimizing the trajectory using the current values
       newtraj(1,i) = traj(1,i-1)-gamma*gx(r,c) ; 
       newtraj(2,i) = traj(2,i-1)-gamma*gy(r,c) ; 
    end
    
    for i=1:length(newtraj(1,:)) 
        newpath_values(i)=obs_cost(floor(newtraj(1,i)),floor(newtraj(2,i)));
    end
    
    % check end conditions
    e = newpath_values - path_values ; 
    mse = 0.5 * sum(e.^2) ; % mean square error
    
    traj = newtraj ; 
    path_values = newpath_values ; 
    iter = iter + 1 ;  
end

% convert my new trajectory to the path which gets plotted in the next
% section
path = newtraj' ; 

%% part (b) 
% 
% gamma = 0.1 ; % provided in question 
% maxiters = 100 ; % max number of iterations 
% eps = 0.1 ; % termination tolerance
% 
% path_values = path_init_values' ; 
% mse = 0.5 * sum(path_values.^2) ; 
% iter = 1 ; 
% 
% % run a while loop until convergence
% while iter < maxiters && mse > eps 
% 
%     % new trajectory
%     newtraj = zeros(2,waypoints);
%     newtraj(1,1)=SX;
%     newtraj(2,1)=SY;
% 
%     for i = 2:waypoints
%        % find the gradient for the current point in the trajectory
%        r = floor(traj(1,i-1)) ; % x position
%        c = floor(traj(2,i-1)) ; % y position
%        if r == 0 
%            r = 1 ; 
%        end 
%        if c == 0 
%            c = 1 ; 
%        end
% 
%        % I am optimizing the trajectory using the current values
%        
%        % smoothness cost
%        e_prev = [traj(1,i-1), traj(2,i-1)] ; 
%        e_now = [traj(1,i), traj(2,i)] ; 
%        smooth = 0.5*(norm(e_now - e_prev))^2 ; 
%        
%        newtraj(1,i) = traj(1,i-1) - gamma * (0.8*gx(r,c) + 4*smooth) ; 
%        newtraj(2,i) = traj(2,i-1) - gamma * (0.8*gy(r,c) + 4*smooth) ; 
%        
%        if newtraj(1,i) < 0 
%            newtraj(1,i) = 0 ;
%        end
%        if newtraj(2,i) < 0
%            newtraj(2,i) = 0 ; 
%        end
%     end
%     
%     for i=1:length(newtraj(1,:)) 
%         r = floor(newtraj(1,i)) ; 
%         if r == 0 
%             r = 1 ; 
%         end 
%         c = floor(newtraj(2,i)) ; 
%         if c == 0 
%             c = 1 ; 
%         end
%         
%         newpath_values(i)=obs_cost(r,c);
%     end
%     
%     % check end conditions
%     e = newpath_values - path_values ; 
%     mse = 0.5 * sum(e.^2) ; % mean square error
%     
%     traj = newtraj ; 
%     path_values = newpath_values ; 
%     iter = iter + 1 ;  
% end
% 
% % convert my new trajectory to the path which gets plotted in the next
% % section
% path = newtraj' ; 

%% part (c) 
% 
% gamma = 0.1 ; % provided in question 
% maxiters = 100 ; % max number of iterations 
% eps = 0.1 ; % termination tolerance
% 
% path_values = path_init_values' ; 
% mse = 0.5 * sum(path_values.^2) ; 
% iter = 1 ; 
% 
% % run a while loop until convergence
% while iter < maxiters && mse > eps 
% 
%     % new trajectory
%     newtraj = zeros(2,waypoints);
%     newtraj(1,1)=SX;
%     newtraj(2,1)=SY;
% 
%     for i = 2:(waypoints-1)
%        % find the gradient for the current point in the trajectory
%        r = floor(traj(1,i-1)) ; % x position
%        c = floor(traj(2,i-1)) ; % y position
%        if r == 0 
%            r = 1 ; 
%        end 
%        if c == 0 
%            c = 1 ; 
%        end
% 
%        % I am optimizing the trajectory using the current values
%        
%        % smoothness cost
%        e_prev = [traj(1,i-1), traj(2,i-1)] ; 
%        e_now = [traj(1,i), traj(2,i)] ; 
%        e_next = [traj(1,i+1), traj(2,i+1)] ; 
%        smooth = 0.5*(norm(e_now - e_prev))^2 + 0.5*(norm(e_next - e_now))^2 ; 
%        
%        newtraj(1,i) = traj(1,i-1) - gamma * (0.8*gx(r,c) + 4*smooth) ; 
%        newtraj(2,i) = traj(2,i-1) - gamma * (0.8*gy(r,c) + 4*smooth) ; 
%        
%        if newtraj(1,i) < 0 
%            newtraj(1,i) = 0 ;
%        end
%        if newtraj(2,i) < 0
%            newtraj(2,i) = 0 ; 
%        end
%     end
%     
%     for i=1:length(newtraj(1,:)) 
%         r = floor(newtraj(1,i)) ; 
%         if r == 0 
%             r = 1 ; 
%         end 
%         c = floor(newtraj(2,i)) ; 
%         if c == 0 
%             c = 1 ; 
%         end
%         
%         newpath_values(i)=obs_cost(r,c);
%     end
%     
%     % check end conditions
%     e = newpath_values - path_values ; 
%     mse = 0.5 * sum(e.^2) ; % mean square error
%     
%     traj = newtraj ; 
%     path_values = newpath_values ; 
%     iter = iter + 1 ;  
% end
% 
% % convert my new trajectory to the path which gets plotted in the next
% % section
% path = newtraj' ; 


%% plot the trajectories
path_values = zeros(tt,1);
for i=1:tt
    r = floor(path(i,1)) ; 
    if r == 0 
        r = 1 ; 
    end 
    c = floor(path(i,2)) ; 
    if c == 0 
        c = 1 ; 
    end
    path_values(i)=obs_cost(r,c);
end
figure(1)
hold on;
plot3(path(:,1),path(:,2),path_values,'.g','MarkerSize',20);

hold off;