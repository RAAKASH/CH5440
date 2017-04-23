function [ G] = gauss( x1,x2,w )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
G = exp(-(x1-x2)'*(x1-x2)/w);

end

