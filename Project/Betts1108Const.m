function [ out ] = Betts1108Const( X )
%Betts8Const Summary of this function goes here
%   Detailed explanation goes here
X = X(:);
if(length(X) ~= 3)
    out = NaN;
    return;
end

out = X(1) + X(2) + X(3) - 1;

end

