function [ out ] = Betts1222Fun( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);
if(length(X) ~= 3)
    out = NaN;
    return;
end

out = - X(1) * X(2) * X(3);

end

