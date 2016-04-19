function [Xopt , Steps] = AugmentedLagrangien(varargin)
%AugmentedLagrangien For Equalities and Inequalities constrant optimization
%  minimize f(x) + ...
% (penalty / 2) * ||h(x)||^2+ lambda' * h(x) + ...
% (penalty / 2) * ||g+(x)||^2+ mu' * g+(x);

% At least 4 variablees
if(nargin < 4)
    error('myfuns:AugmentedLagrangien:WrongInput', ...
        'Should at least give the functon to optimize f, and equality constaints h, and initial conditions X0, L0');
end

if(isa(varargin{3},'function_handle'))
    inequality = 1;
    f = varargin{1};
    h = varargin{2};
    g = varargin{3};
    X0 = varargin{4};
    L0 = varargin{5};
    M0 = varargin{6};
else
    inequality = 0;
    f = varargin{1};
    h = varargin{2};
    X0 = varargin{3};
    L0 = varargin{4};
end

if(isa(h,'function_handle'))
    equality = 1;
else
    equality = 0;
end

% default option
N = 100;
method = 'bfgs';
psi = 1.5;
r = 1;
epsX = 1e-4;
epsC = 1e-8;
verbose = 0;
df = [];
dh = [];
dg = [];
dv = [];

% options
for nVar = (4 + 2 * inequality + 1):2:nargin
    switch(lower(varargin{nVar}))
        case 'maxiter'
            N = varargin{nVar+1};
        case 'method'
            method = lower(varargin{nVar+1});
        case 'penaltymulti'
            psi = varargin{nVar+1};
        case 'penalty'
            r = varargin{nVar+1};
        case 'epsstate'
            epsX = varargin{nVar+1};
        case 'epscost'
            epsC = varargin{nVar+1};            
        case 'fdiff'
            df = varargin{nVar+1};
        case 'hdiff'
            dh = varargin{nVar+1};
        case 'gdiff'
            dg = varargin{nVar+1};
        case 'verbose'
            switch(lower(varargin{nVar+1}))
                case {'testing', 'test', 't', 3}
                    verbose = 3;
                case {'iterative', 'iter', 'i', 2}
                    verbose = 2;
                case {'notify', 'yes', 'y', 1}
                    verbose = 1;
                otherwise
                    verbose = 0;
            end
        otherwise
            error('myfuns:AugmentedLagrangien:WrongInput', ...
                'Unkown option %s', varargin{nVar});
    end
end

X0 = X0(:);
L0 = L0(:);

Steps.X = NaN*ones(length(X0),N);
Steps.f = NaN*ones(1, N);
if(equality)
    Steps.l = NaN*ones(length(L0),N);
    Steps.h = NaN*ones(length(L0),N);
end
if(inequality)
    Steps.m = NaN*ones(length(M0),N);
    Steps.g = NaN*ones(length(M0),N);
end
Steps.r = NaN*ones(1,N);
Steps.psi = NaN*ones(1,N);
Steps.Time = 0;

n = 1;
Steps.X(:,n) = X0;
Steps.f(n) = f(Steps.X(:,n));
if(equality)
    Steps.l(:,n) = L0;
    Steps.h(:,n) = h(Steps.X(:,n));
end
if(inequality)
    Steps.m(:,n) = M0;
    Steps.g(:,n) = g(Steps.X(:,n));
end
Steps.r(n) = r;
Steps.psi(n) = psi;


LineVerbose{n} = sprintf(               '+---------------------------------------------------------------- Augmented Lagrangien, Max iterations %5d: -----------------------------------------------------+\\n', N);
LineVerbose{n} = strcat(LineVerbose{n}, '$ Iter |                  State                   |   Function    |           Equality Constraints           |           Inequality Constraints         |  Penalty $\n');
LineVerbose{n} = strcat(LineVerbose{n}, '+-------------------------------------------------+---------------+------------------------------------------+------------------------------------------+----------+\n');
if(inequality)
    if(equality)
        LineVerbose{n} = strcat(LineVerbose{n}, sprintf('$%06d| %40s | %+8.6e | %40s | %40s | %8.2f $ \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.h(:,n)), Vector2String(Steps.g(:,n)), Steps.r(:,n)));
    else
        LineVerbose{n} = strcat(LineVerbose{n}, sprintf('$%06d| %40s | %+8.6e | ---------------------------------------- | %40s | %8.2f $ \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.g(:,n)), Steps.r(:,n)));
    end
else
    LineVerbose{n} = strcat(LineVerbose{n}, sprintf('$%06d| %40s | %+8.6e | %40s | ---------------------------------------- | %8.2f $ \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.h(:,n)), Steps.r(:,n)));
end

if(verbose > 0)
    fprintf(LineVerbose{n});
end

while n < N
    
    rn = Steps.r(n);
    if(equality)
        ln = Steps.l(:,n);
    end
    if(inequality)
        mn = Steps.m(:,n);
    end
    
    if(inequality)
        gp = @(x)(max(g(x),-mn/rn));
        if(equality)
            v = @(x)(f(x) + ...
                rn * (h(x)' * h(x)) / 2 + ...
                ln' * h(x) + ...
                rn * (gp(x)' * gp(x)) / 2 + ...
                mn' * gp(x));
            if(~isempty(df) && ~isempty(dh) && ~isempty(dg))
                dgp = @(x)( ((ones(length(X0), 1)) * (g(x) > -mn/rn)').*dg(x) );
                
                dv = @(x)(df(x) + ...
                    dh(x) * (ln + rn*h(x)) + ...
                    dgp(x) * (mn + rn*g(x)));
            end
        else
            v = @(x)(f(x) + ...
                rn * (gp(x)' * gp(x)) / 2 + ...
                mn' * gp(x));
            if(~isempty(df) && ~isempty(dg))
                dgp = @(x)( ((ones(length(X0), 1)) * (g(x) > -mn/rn)').*dg(x) );
                
                dv = @(x)(df(x) + ...
                    dgp(x) * (mn + rn*g(x)));
            end
        end
    else
        v = @(x)(f(x) + ...
            rn * (h(x)' * h(x)) / 2+ ...
            ln' * h(x));
        if(~isempty(df) && ~isempty(dh))
            dv = @(x)(df(x) + ...
                dh(x) * (ln + rn*h(x)));
        end
    end
    
    [Steps.X(:,n+1), Steps.secant{n}] = Secant(v, dv, Steps.X(:,n), 'MaxIter', N*1e3, 'Method', method, 'epsState', epsX, 'epsCost', epsC, 'verbose', verbose - 1);
    Steps.Time = Steps.Time + Steps.secant{n}.Time;
    
    if(isnan(Steps.X(:,n+1)))
        error('myfuns:AugmentedLagrangien:NaN', ...
            'Error in computations');
    end
    
    Steps.psi(n + 1) = Steps.psi(n);
    
    % Stop Conditions
    if(equality && max(abs(h(Steps.X(:,n+1)))) < epsC)
        if(~inequality)
            break;
        elseif(inequality && sum(g(Steps.X(:,n+1)) > epsC) == 0)
            break;
        else % if inequalities are not satisfied yet
            % continue
        end
    end
    
    if(max(abs(Steps.X(:,n+1) - Steps.X(:,n))) < epsX)
        if(~inequality)
            break;
        elseif(inequality && sum(g(Steps.X(:,n+1)) > epsC) == 0)
            break;
        else % if inequalities are not satisfied yet
            % continue
        end
    end
    
    Steps.r(n+1) = Steps.psi(n + 1) * rn;
    if(equality)
        Steps.l(:,n+1) = ln + rn*h(Steps.X(:,n+1));
    end
    if(inequality)
        Steps.m(:,n+1) = max(0, mn + rn*g(Steps.X(:,n+1)));
    end
    
    n = n+1;
    
    Steps.f(n) = f(Steps.X(:,n));
    if(equality)
        Steps.h(:,n) = h(Steps.X(:,n));
    end
    if(inequality)
        Steps.g(:,n) = g(Steps.X(:,n));
    end
    
    if(inequality)
        if(equality)
            LineVerbose{n} = sprintf('$%06d| %40s | %+8.6e | %40s | %40s | %8.2f $ \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.h(:,n)), Vector2String(Steps.g(:,n)), Steps.r(:,n));
        else
        LineVerbose{n} = sprintf('$%06d| %40s | %+8.6e | ---------------------------------------- | %40s | %8.2f $ \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.g(:,n)), Steps.r(:,n));
        end
    else
        LineVerbose{n} = sprintf('$%06d| %40s | %+8.6e | %40s | ---------------------------------------- | %8.2f $ \\n', n, Vector2String(Steps.X(:,n)), Steps.f(n), Vector2String(Steps.h(:,n)), Steps.r(:,n));
    end
    
    if(verbose > 0)
        fprintf(LineVerbose{n});
    end
    
end

if(verbose > 0)
    LineVerbose{n+1} = '+-------------------------------------------------+---------------+------------------------------------------+------------------------------------------+----------+\n';
    fprintf(LineVerbose{n+1});
end

Xopt = Steps.X(:,n);
if(sum(sum(isnan(Steps.X))))
    Steps.X(isnan(Steps.X)) = [];
    Steps.X = reshape(Steps.X, length(X0), length(Steps.X)/length(X0));
    Steps.f(isnan(Steps.f)) = [];
    Steps.r(isnan(Steps.r)) = [];
    if(equality)
        Steps.h(isnan(Steps.h)) = [];
        Steps.h = reshape(Steps.h, length(L0), length(Steps.h)/length(L0));
        Steps.l(isnan(Steps.l)) = [];
        Steps.l = reshape(Steps.l, length(L0), length(Steps.l)/length(L0));
    end
    if(inequality)
        Steps.g(isnan(Steps.g)) = [];
        Steps.g = reshape(Steps.g, length(M0), length(Steps.g)/length(M0));
        Steps.m(isnan(Steps.m)) = [];
        Steps.m = reshape(Steps.m, length(M0), length(Steps.m)/length(M0));
    end
end
Steps.NumOfIterations = n;

Steps.Verbose = sprintf(cell2mat(LineVerbose));

end

