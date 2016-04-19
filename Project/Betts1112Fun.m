function [ out ] = Betts1112Fun( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);
if(length(X) ~= 2)
    out = NaN;
    return;
end

out = (1 - X(1))^2;

end

