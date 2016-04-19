function [ out ] = Betts1208Const( X )
%Betts8Const Summary of this function goes here
%   Detailed explanation goes here
X = X(:);
if(length(X) ~= 2)
    out = NaN;
    return;
end

out(1) = -X(1) - X(2)^2;
out(2) = -X(1)^2 - X(2);
out(3) = X(1) - 0.5;
out(4) = -X(1) - 0.5;
out(5) = -X(1)^2 - X(2)^2 + 1;

out = out(:);

end

