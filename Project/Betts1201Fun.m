function [ out ] = Betts1201Fun( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);
if(length(X) ~= 5)
    out = NaN;
    return;
end

out = 2 - X(1) * X(2) * X(3) * X(4) * X(5) / 120;

end

