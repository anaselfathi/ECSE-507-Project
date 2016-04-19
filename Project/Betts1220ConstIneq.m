function [ out ] = Betts1220ConstIneq( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);

if(length(X) ~= 16)
    out = NaN;
    return;
end

out = [-X; X - 5];

out = out(:);

end

