function [con] = Betts1102Const(x)
con = x(1)*(1+x(2)^2)+(x(3)^4)-4-(3*sqrt(2));
end