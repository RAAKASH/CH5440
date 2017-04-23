function [ Psat] = p( T )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
A = 14.0568;
B = 2825.42;
C = 230.44;
Psat = exp(A - B./(T+C));
end

