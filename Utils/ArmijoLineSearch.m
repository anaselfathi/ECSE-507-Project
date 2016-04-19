function w = ArmijoLineSearch(v, dvx, x, s, gamma, mu, N)

% Check number of inputs.
if nargin > 7
    error('myfuns:armijoLineSearch:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% Fill in unset optional values.
switch nargin
    case 4
        gamma = 1.5;
        mu = 0.8;
        N = 100;
    case 5
        mu = 0.8;
        N = 100;
    case 6
        N = 100;
end

x = x(:);
s = s(:);
dvx = dvx(:);
vx = v(x);

n = length(x);

if n ~= length(s)
    error('myfuns:armijoLineSearch:TooManyInputs', ...
        's and x should have the same length');
end

if n ~= length(dvx)
    error('myfuns:armijoLineSearch:TooManyInputs', ...
        'dvx and x should have the same length');
end

p = 1;
q = 1;

% Find a w on right of w^
w = 1;
while(p < N)
    w = w * gamma;
    if(v(x + w * s) > vx + 0.5 * w * (dvx'*s))
        break;
    end
    p = p + 1;
end

% Find a w on left of w^
while(q < N)
    w = w * mu;
    if(v(x + w * s) < vx + 0.5 * w * (dvx'*s))
        break;
    end
    q = q + 1;
end

end