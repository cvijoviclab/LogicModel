function edgepoints = rectangleMaker(c,a,b)

% edgepoints will contain the four edge points of the rectangle with 
% length a and width b defined by the center c
% the points are order in a counterclockwise fashion and the first point
% was added at the end again so that the plot function actually plots the
% whole rectangle
% size of edgepoints: 5x2

    edge1 = c + [a/2;b/2];
    edge2 = c + [-a/2;+b/2];
    edge3 = c + [-a/2;-b/2];
    edge4 = c + [a/2;-b/2];

    edgepoints = [edge1 edge2 edge3 edge4 edge1]';

end