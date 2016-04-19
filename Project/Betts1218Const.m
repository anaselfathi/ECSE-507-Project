function [ out ] = Betts1218Const( X )
%BETTS6FUN Summary of this function goes here
%   Detailed explanation goes here

X = X(:);

if(length(X) ~= 5)
    out = NaN;
    return;
end

out(1) = 85.334407 + 0.0056858 * X(2) * X(5) + 0.0006262*X(1)*X(4)- 0.0022053*X(3)*X(5);
out(2) = 92 - out(1);
out(3) = 80.51249 + 0.0071317*X(2)*X(5) + 0.0029955*X(1)*X(2) + 0.0021813*X(3)^2- 90;
out(4) = 20 - out(3);
out(5) = 9.300961 + 0.0047026*X(3)*X(5) + 0.0012547*X(1)*X(3) + 0.0019085*X(3)*X(4)-20;
out(6) = 5 - out(5);
out(7) = X(1) - 78;
out(8) = 102 - X(1);
out(9) = X(2) - 33;
out(10) = 45 - X(2);
out(11) = X(3)- 27;
out(12) = 45 - X(3);
out(13) = X(4) - 27;
out(14) = 45 - X(4);
out(15) = X(5) - 27;
out(16) = 45 - X(5);

out = -out(:);

end

