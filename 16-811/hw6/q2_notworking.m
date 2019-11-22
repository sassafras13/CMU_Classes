% 16-811 Fall 2019
% Emma Benjaminson
% Assignment 6
% Problem 2
% References: 
% [1] https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
% [2] http://www.cs.kent.edu/~dragan/ST-Spring2016/visibility%20graphs.pdf
% [3] https://medium.com/basecs/finding-the-shortest-path-with-a-little-help-from-dijkstra-613149fbdc8e
% [4] https://brilliant.org/wiki/dijkstras-short-path-finder/

%%
clc ; clear all ; close all ;

%% Main

%%%%%%%%%%%%%%%%%%%%%%%%
% START AND END POINTS %
%%%%%%%%%%%%%%%%%%%%%%%%
% take in the start point and the goal point
start = [1,1] ; 
goal = [120,120] ; 

figure(1)
plot(start(1), start(2), 'or') 
hold on 
plot(goal(1), goal(2), 'or') 
hold on 

%%%%%%%%%%%%%
% OBSTACLES %
%%%%%%%%%%%%%

% generate a list of obstacles which are convex polygons
No = 2 ; % number of obstacles
maxV = 10 ; % maximum number of vertices permitted in one obstacle
maxL = 50 ; % maximum length of an obstacle

obstacles = {} ; % initialize an empty cell array that can contain all obstacle data

for i = 1:No
    Nv = 1 ; % need to set a condition to force number of vertices to be at least 3
    
    while Nv < 3
        Nv = round(maxV*rand(1)) ; % number of vertices on the obstacle
        offset = round(100*rand(1)) ; % offset of location of convex obstacle
    end
    
    finalx = convexObstacle(Nv,maxL,offset) ; % produce a list of points 
    pgon = polyshape(finalx(:,1),finalx(:,2),'Simplify',false) ; 
    obstacles(i,:) = {finalx} ; % add the list of points to the obstacles cell array
    
%     color = [rand(1), rand(1), rand(1)] ; % set color at random for this obstacle
    
    figure(1)
    plot(pgon)
    hold on
%     plot(finalx(:,1), finalx(:,2), 'o', 'color', color) 
%     hold on
%     for i = 1:(length(finalx)-1)
%         plot(finalx(i:i+1, 1), finalx(i:i+1, 2), '-', 'color', color)
%     end
%     hold on
%     plot([finalx(end,1),finalx(1,1)], [finalx(end,2),finalx(1,2)],'-', 'color', color)
%     hold on

end

% obstacles = [   44.3019   52.0579 ;
%    36.9086   46.3763 ;
%    53.0475   40.3065 ] ; 
% No = 1 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF START/END IS INSIDE OBSTACLE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:No
    obs = obstacles{i} ; 
    
    alertStart = isContained(start,obs) ; 
    alertGoal = isContained(goal,obs) ; 
    
    if alertStart == 1 
        error('Start point is contained inside an obstacle, ending script.') ;
    end
    
    if alertGoal == 1
        error('Goal point is contained inside an obstacle, ending script.') ; 
    end
end

%%%%%%%%%%%%%%%%%%%%%
% VISIBILITY GRAPHS %
%%%%%%%%%%%%%%%%%%%%%

% generate visGraph for start point
visGraphStart = visGraphGen(start, obstacles) ; 

n = size(visGraphStart,1) ; 

% for i = 1:n
%     figure(1)
%     plot([visGraphStart(i,1),visGraphStart(i,3)],[visGraphStart(i,2),visGraphStart(i,4)],'-k')
%     hold on
% end

% generate visGraph between each obstacle
visGraphObs = [] ; 

% run through each obstacle
for i = 1:No
    obstacle = obstacles{i} ; 
%     obstacle = obstacles ; 
    vertices = size(obstacle,1) ; 
    
    % for each vertex in each obstacle, get the visibility graph for that
    % vertex and add to the list of visibility graph for all obstacles
    for j = 1:vertices
        vG = visGraphGenObs(obstacle(j,:),obstacles) ; 
        visGraphObs = [visGraphObs ; vG] ; 
    end
end

n = size(visGraphObs,1) ; 

% for i = 1:n
%     figure(1)
%     plot([visGraphObs(i,1),visGraphObs(i,3)],[visGraphObs(i,2),visGraphObs(i,4)],'-k')
%     hold on
% end

% generate visGraph for goal point
visGraphEnd = visGraphGen(goal, obstacles) ; 

temp = visGraphEnd ; 
visGraphEnd(:,1:2) = temp(:,3:4) ; 
visGraphEnd(:,3:4) = temp(:,1:2) ; 

n = size(visGraphEnd,1) ; 

% for i = 1:n
%     figure(1)
%     plot([visGraphEnd(i,1),visGraphEnd(i,3)],[visGraphEnd(i,2),visGraphEnd(i,4)],'-k')
%     hold on
% end


%%%%%%%%%%%%%%%%%%%%%%%%
% DIJKSTRA'S ALGORITHM %
%%%%%%%%%%%%%%%%%%%%%%%%

% put all the visibility graphs into one variable
Graph = [visGraphStart ; visGraphObs; visGraphEnd ] ; 

% get the number of nodes in the graph
% build a list of all the nodes
Nn = 2 ; 
Nodes = [start] ; 

for i = 1:No
    n = length(obstacles{i}) ;
%     n = 3 ; 
    Nn = Nn + n ;
    Nodes = [Nodes ; obstacles{i}] ;
%     Nodes = [Nodes ; obstacles] ; 
end

Nodes = [Nodes ; goal] ; 

% figure(1)
% for i = 1:length(Graph)
%     plot([Graph(i,1),Graph(i,3)],[Graph(i,2),Graph(i,4)],'-k') ; 
%     hold on
% end
% 
% for i = 1:length(Nodes)
%     plot(Nodes(i,1),Nodes(i,2),'or') 
%     hold on
% end

[S, Q] = dijkstra(Nodes, Graph) ; 

% use S to calculate the shortest distance from goal to start
waypoints = shortestPath(S,start,goal) ; 

for i = 1:(length(waypoints)-1)
    figure(1)
    plot(waypoints(i:i+1,1), waypoints(i:i+1,2), '-r','LineWidth',3)
    hold on
end

figure(1)
title('HW3 Q2') ; xlabel('X') ; ylabel('Y') ; 

%% Extra Code

% % test case for intersection function
% L1 = [1, 1 ; 5, 4] ; 
% L2 = [2, 3 ; 5, 1] ; 
% L3 = [1, 4 ; 5, 7] ; 
% result1 = intersection(L1, L2) ;
% result2 = intersection(L2, L3) ;

%%% TEST CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test1 = [36.3756   20.3943 ;
%    30.2903   32.0053 ;
%    18.1838   24.5737 ;
%    17.4860   14.2309 ;
%    24.5264   12.0560 ] ;
% 
%     color = [rand(1), rand(1), rand(1)] ; % set color at random for this obstacle
%     
% figure(1)
% 
% plot(test1(:,1), test1(:,2), 'o', 'color', color) 
% hold on
% for i = 1:(length(test1)-1)
%     plot(test1(i:i+1, 1), test1(i:i+1, 2), '-', 'color', color)
% end
% hold on
% plot([test1(end,1),test1(1,1)], [test1(end,2),test1(1,2)],'-', 'color', color)
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions

function finalx = convexObstacle(Nv,maxL,offset)
    
    % run a function to generate a set of n random points in 2D space
    x = randpoints(Nv,maxL,offset) ; 

    % look for the point that has the smallest y value
    % if there are multiple such points, choose the one with the largest x
    % value
    [originy, ii] = min(x(:,2)) ; % ii is the index of the minimum value

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
    stack = [] ;
    stack = [stack ; x(1:2,:)] ;

    for i = 3:length(x) 
        while length(stack) > 1 && isconvex(stack(end-1,1:2),stack(end,1:2),x(i,1:2)) == 0 
            stack(end,:) = [] ; 
        end
        stack = [stack ; x(i,:)] ; 
    end

    finalx = stack(:,1:2) ; 

end

function x = randpoints(Nv,maxL,offset)
    x = maxL*rand([Nv,2]) + offset ; % scale the random numbers so they are between 0 and 100
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

function result = intersection(L1, L2)
    
    eps = 0.2 ; 

    % break lines into [(x1, y1), (x2, y2)] and [(x3, y3) , (x4, y4)] coordinates
    x1 = L1(1) ; 
    y1 = L1(2) ; 
    x2 = L1(3) ; 
    y2 = L1(4) ; 
    
    x3 = L2(1) ; 
    y3 = L2(2) ; 
    x4 = L2(3) ; 
    y4 = L2(4) ; 
    
    % calculate the slopes and intercepts of the two lines 
    m1 = (y2 - y1) / (x2 - x1) ; 
    m2 = (y4 - y3) / (x4 - x3) ; 
    
    b1 = y1 - m1*x1 ; 
    b2 = y3 - m2*x3 ; 
    
    % check first that there exists some interval in x where an
    % intersection could exist
    if max(x1,x2) < min(x3,x4)
        result = false ; 
        return % if this result is false, end the function
    end
    
    % check if the slopes of the lines are equal, then the lines are
    % parallel and there is not point of intersection
    % set the difference between the two slopes to be very small, but
    % non-zero, in the event that they are not exactly the same value
    if abs(m1 - m2) < 0.001
        result = false ; 
        return
    end
    
    % the potential point of intersection (x,y) must lie inside the
    % intervals I1 and I2
    Ix1 = [min(x1,x2), max(x1,x2)] ; 
    Ix2 = [min(x3,x4), max(x3,x4)] ; 
    
    Iy1 = [min(y1,y2), max(y1,y2)] ; 
    Iy2 = [min(y3,y4), max(y3,y4)] ; 
    
    % and the potential point of intersection really needs to fall in the
    % intersection of these two intervals
    Ixa = [max(Ix1(1),Ix2(1)), min(Ix1(2),Ix2(2))] ; 
    Iya = [max(Iy1(1),Iy2(1)), min(Iy1(2),Iy2(2))] ; 
    
    % if there is a point of intersection its x value will be
    xa = (b2 - b1) / (m1 - m2) ; 
    xa = xa + 0.1*rand(1) ; 
    
    % and its y value will be
    ya = m1*xa + b1 ; 
    ya = ya + 0.1*rand(1) ; 
    
    % now check that xa lies within the interval Ia
    if ((xa - eps) < Ixa(1) || (xa + eps) > Ixa(2)) && ((ya - eps) < Iya(1) || (ya + eps) > Iya(2))
        result = false ; 
        return
    else
        result = true ; 
    end
end

function result = intersectionObs(L1, L2)
    
    eps = 0.2 ; 

    % break lines into [(x1, y1), (x2, y2)] and [(x3, y3) , (x4, y4)] coordinates
    x1 = L1(1) ; 
    y1 = L1(2) ; 
    x2 = L1(3) ; 
    y2 = L1(4) ; 
    
    x3 = L2(1) ; 
    y3 = L2(2) ; 
    x4 = L2(3) ; 
    y4 = L2(4) ; 
    
    % calculate the slopes and intercepts of the two lines 
    m1 = (y2 - y1) / (x2 - x1) ; 
    m2 = (y4 - y3) / (x4 - x3) ; 
    
    b1 = y1 - m1*x1 ; 
    b2 = y3 - m2*x3 ; 
    
    % check first that there exists some interval in x where an
    % intersection could exist
    if max(x1,x2) < min(x3,x4)
        result = false ; 
        return % if this result is false, end the function
    end
    
    % check if the slopes of the lines are equal, then the lines are
    % parallel and there is not point of intersection
    % set the difference between the two slopes to be very small, but
    % non-zero, in the event that they are not exactly the same value
    if abs(m1 - m2) < 0.001
        result = false ; 
        return
    end
    
    % the potential point of intersection (x,y) must lie inside the
    % intervals I1 and I2
    Ix1 = [min(x1,x2), max(x1,x2)] ; 
    Ix2 = [min(x3,x4), max(x3,x4)] ; 
    
    Iy1 = [min(y1,y2), max(y1,y2)] ; 
    Iy2 = [min(y3,y4), max(y3,y4)] ; 
    
    % and the potential point of intersection really needs to fall in the
    % intersection of these two intervals
    Ixa = [max(Ix1(1),Ix2(1)), min(Ix1(2),Ix2(2))] ; 
    Iya = [max(Iy1(1),Iy2(1)), min(Iy1(2),Iy2(2))] ; 
    
    % if there is a point of intersection its x value will be
    xa = (b2 - b1) / (m1 - m2) ; 
    xa = xa + 1*rand(1) ; 
    
    % and its y value will be
    ya = m1*xa + b1 ; 
    ya = ya + 1*rand(1) ; 
    
    % now check that xa lies within the interval Ia
    if ((xa - eps) < Ixa(1) || (xa + eps) > Ixa(2)) && ((ya - eps) < Iya(1) || (ya + eps) > Iya(2))
        result = false ; 
        return
    else
        result = true ; 
    end
end

function edges = findedges(obstacle) 
    n = length(obstacle) ; % number of vertices (and also edges)
        
    edges = zeros(n,4) ; % each edge is described by 4 points [x1 y1 x2 y2]
    
    for i = 1:(n-1)
        edges(i,1:2) = obstacle(i,:) ; % start point of edge
        edges(i,3:4) = obstacle(i+1,:) ; % end point of edge
    end
    
    % the last edge runs from the last vertex to the first one, add to list
    edges(n,1:2) = obstacle(n,:) ; 
    edges(n,3:4) = obstacle(1,:) ; 
    
end

function visGraph = visGraphGenObs(point, obstacles) 
    eps = 1 ; 
    No = size(obstacles,1) ; 
    obstacle = obstacles ; 
    visGraph = [] ; 

    % find all the edges of all the obstacles
    edges = [] ; 
%     edges = [edges; findedges(obstacle)] ;

    for i = 1:No
        % pull out the data for this obstacle from the cell
        obstacle = obstacles{i} ; 

        % find the edges for this obstacle and add to array
        edges = [edges; findedges(obstacle)] ; 
    end

    Ne = size(edges,1) ; 

    for i = 1:No
        obstacle = obstacles{i} ; 
        Nv = size(obstacle,1) ; % number of vertices in the obstacle

        % for every vertex in the obstacle
        for j = 1:Nv

            % join the start point to the jth vertex
            L1 = [point, obstacle(j,:)] ;

            results = zeros(Ne,1) ; 

            % check every edge for intersection with L1
            for k = 1:Ne

                L2 = edges(k,:) ; 

                % check if the line segments intersect the edges of the obstacle
                results(k) = intersectionObs(L1, L2) ; 

            end

            % if L1 does not intersect any edge, add L1 to visGraph
            if sum(results) == 0 
                visGraph = [visGraph ; L1] ; 
            end
        end
    end
end

function visGraph = visGraphGen(point, obstacles) 
    eps = 1 ; 
    No = size(obstacles,1) ; 
    obstacle = obstacles ; 
    visGraph = [] ; 

    % find all the edges of all the obstacles
    edges = [] ; 
%     edges = [edges; findedges(obstacle)] ;

    for i = 1:No
        % pull out the data for this obstacle from the cell
        obstacle = obstacles{i} ; 

        % find the edges for this obstacle and add to array
        edges = [edges; findedges(obstacle)] ; 
    end

    Ne = size(edges,1) ; 

    for i = 1:No
        obstacle = obstacles{i} ; 
        Nv = size(obstacle,1) ; % number of vertices in the obstacle

        % for every vertex in the obstacle
        for j = 1:Nv

            % join the start point to the jth vertex
            L1 = [point, obstacle(j,:)] ; 

%             figure(2)
%             plot([L1(1),L1(3)],[L1(2),L1(4)],'-r')
%             hold on

            results = zeros(Ne,1) ; 

            % check every edge for intersection with L1
            for k = 1:Ne

                L2 = edges(k,:) ; 

                % check if the line segments intersect the edges of the obstacle
                results(k) = intersection(L1, L2) ; 

            end

            % if L1 does not intersect any edge, add L1 to visGraph
            if sum(results) == 0 
                visGraph = [visGraph ; L1] ; 
            end
        end
    end
end

function dist = distanceNodes(L)
    x1 = L(1) ; 
    y1 = L(2) ; 
    x2 = L(3) ;
    y2 = L(4) ; 
    
    dx = abs(x2 - x1) ; 
    dy = abs(y2 - y1) ; 
    
    dist = sqrt((dx^2) + (dy^2)) ; 
end

function [S, Q] = dijkstra(Nodes, Graph)
    % create a variable that contains the distances to each node from the
    % source
    dist = zeros(length(Nodes),1) ; 

    % set the distance to the source node as 0 and set the distance to all other nodes as infinity
    for i = 2:length(Nodes)
        dist(i) = Inf ; 
    end

    % add each node to a list Q
    Q = [Nodes, dist, zeros(length(Nodes),2)] ; 

    % create list S of visited nodes 
    S = [] ; 

    % while Q is not empty
    while isempty(Q) == 0

        % find the node v with the minimum distance and remove
        % that node from the list Q
        [~,vind] = min(Q(:,3)) ; 
        v = Q(vind,1:2) ; 

        % for every neighbor of the node v that is still in Q
        u = [] ; 

        for i = 1:length(Graph)
            if Graph(i,1:2) == v 
                u = [u ; Graph(i,3:4)] ; 
            end
        end

        % calculate the distance from v to the neighbor, u as the distance of v from the source node plus
        % the distance between v and u
        for i = 1:size(u,1)
            if ismember(u(i,:),S) == 0
                L = [v , u(i,:)] ; 
                distN = distanceNodes(L) ;
                newDist = distN + Q(vind,3) ; 

                % if this new calculated ddistance is shorter than the assigned distance to
                % u, then update the distance to u with this new calculation
                index = find(Q(:,1:2) == u(i,:),1) ; 
                if newDist < Q(index,3) 
                    Q(index,3) = newDist ; 
                    Q(index,4:5) = v ; 
                end
            end
        end

        if ismember(v(1),S) == 0 
            S = [S ; Q(vind,:)] ; 
            Q(vind,:) = [] ; 
        end
    end
end

function waypoints = shortestPath(S,start,goal)
    cost = 0 ; 
    arrival = [0,0] ; 
    waypoints = [] ; 

    while arrival ~= start
        ii = 1 ;
        while S(ii,1:2) ~= goal
            ii = ii + 1 ; 
        end
        cost = cost + S(ii,3) ; 
        waypoints = [waypoints ; goal] ; 
        arrival = S(ii,1:2) ; 
        goal = S(ii,4:5) ; 
    end
end

function alert = isContained(point,obs)

    % draw a ray extending from the point to the right along the x axis
    L1 = [point, point(1) + 100, point(2)] ; 
    
    % find the edges of the obstacle 
    edges = findedges(obs) ; 
    
    results = zeros(size(edges,1),1) ; 
    
    % check if the ray intersects any of the edges
    for i = 1:size(edges,1)
        L2 = edges(i,:) ; 
        results(i) = intersection(L1,L2) ; 
    end
    
    % count the number of intersections
    num = sum(results) ; 
    
    % return alert == true if the # intersections is odd
    if mod(num,2) == 0 % even
        alert = false ; 
    elseif mod(num,2) == 1 % odd
        alert = true ; 
    end
    
end