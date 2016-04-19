function [ out ] = Betts1108Fun( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);
if(length(X) ~= 3)
    out = NaN;
    return;
end

if(sum(X > 1) ~= 0)
    out = inf;
    return;
end

if(sum(X < 0) ~= 0)
    out = inf;
    return;
end

g = 32.174;
W = 0.03;
a1 = 0.09;
a2 = 0.07;
a3 = 0.13;
I1 = 255;
I2 = 280;
I3 = 290;

out = g * I1 * log((X(1) + X(2) + X(3) + W)/(X(1) * a1 + X(2) + X(3) + W)) + ...
    g * I2 * log((X(2) + X(3) + W)/(X(2) * a2 + X(3) + W)) + ...
    g * I3 * log((X(3) + W)/(X(3) * a3 + W));

out = -out;

end

