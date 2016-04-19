function [ out ] = Betts1222Const( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);

if(length(X) ~= 3)
    out = NaN;
    return;
end

out(1:3) = X;
out(4:6) = 42 - X;
out(7) = X(1) + 2*X(2) + 2*X(3); 
out(8) = 72 - out(7);

out = -out(:);

end

