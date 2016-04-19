function [ out ] = Betts1227FunDiff( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);
if(length(X) ~= 5)
    out = NaN;
    return;
end

a(1) = -24345;
a(2) = -8720288.849;
a(3) = 150512.5253;
a(4) = -156.6950325;
a(5) = 476470.3222;
a(6) = 729482.8271;

out = [a(2) + a(3) * X(2) + a(4) * X(3) + a(5) * X(4) + a(6) * X(5);...
    a(3) * X(1);...
    a(4) * X(1);...
    a(5) * X(1);...
    a(6) * X(1)];

out = -out;

end

