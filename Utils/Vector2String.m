function [ str ] = Vector2String( X )
%VECTOR2STRING Summary of this function goes here
%   Detailed explanation goes here

str = '[';
for i = 1:length(X)
    str = strcat(str,sprintf('%+1.2e, ', X(i)));
end
str = str(1:end-1);
str = strcat(str,']');

end

