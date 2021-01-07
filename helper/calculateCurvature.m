function [curvature] = calculateCurvature(x,y)
%CALCULATECURVATURE Summary of this function goes here
%   Detailed explanation goes here
dx  = gradient(x);
ddx = gradient(dx);
dy  = gradient(y);
ddy = gradient(dy);
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom .* denom .* denom;
curvature = num ./ denom;
%curvature(denom < 0) = NaN;
end

