function [x, y] = ellipseMaker(c, a, b)

% x and y will give 200 points that define the ellipse with major axis a
% and minor axis b
% size of x and y: 200x1

    theta = linspace(0, 2*pi, 200)';
    x = (c(1).*ones(length(theta),1)) + a.*cos(theta);
    y = (c(2).*ones(length(theta),1)) + b.*sin(theta);

end