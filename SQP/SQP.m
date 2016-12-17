% ============================================================================
% This MATLAB function represents non-convex optimization problem solver which
% implements the Sequential quadratic programming method.
% 
% 
%       min  f(x)
%       s.t. h(x) <= 0
%  
% where f(x, y) and also h(x, y) could be convex or non-convex symbolic 
% two-dimensional functions of the variables x and y.
% If we consider more constraints, the h(x,y) will be matrix
% of symbolic expresions which represent left side of constraints.
% 
%       x = SQP(fx, hx, x0, kmax, eps)
%
% Input arguments:
%   
%    fx    - (cost) function with symbolic variables  x, y
%    hx    - matrix (column vector) of symbolic expresions which represent 
%            left side of constraints
%    x0    - Iniatial point (in this phase FEASIBLE!!!)
%    kmax  - maximal no. of iterations for SQP
%    eps   - precision of convergence check
% ============================================================================
function convergence = SQP(fx, hx, x0, kmax, eps)

nc = size(hx,1);    % number of constraints
nx = size(x0,1);    % number of optimization variables

options = sdpsettings('verbose', 0, 'solver', 'quadprog');
convergence = x0;

for k=1:kmax
    grad = gradient(fx);
    hess = hessian(fx);
    
    % ----- approximation variables -----
    x = x0(1);
    y = x0(2);
    
    fx0 = double(subs(fx));
    grad0 = double(subs(grad));
    hess0 = double(subs(hess));
    
    hx0 = double(subs(hx));
    % ---- check if hessian is positive semi-definite ----
    lambda = eig(hess0);
    nl = size(lambda,1);
    positive = 1;
    for i=1:nl
        if (lambda(i) < 0), positive = 0; end
    end                  % the hessian is not positive semi-definite when at least 1 eigenvalue is negative
    if ( positive ~= 1)     
        q = size(hess0);
        hess0 = zeros(q);
    end
    % ---- approximated optimization problem ----
    xx = sdpvar(nx, 1);
    
    cost = fx0 + grad0'*(xx - x0) + (1/2).*(xx - x0)'*hess0*(xx - x0);
    cst = [];
    for i=1:nc
        grad_hx = gradient(hx(i));
        grad_hx0 = double(subs(grad_hx));
        cst = cst + [ hx0(i) + grad_hx0'*(xx - x0) <= 0];
    end
    optimize(cst, cost, options);
    
    % ---- check convergence --------------------
    x_opt = value(xx);
    if ( abs(x_opt - x0) <= eps )
        break
    else
        x0 = x_opt;
        convergence = [convergence , x0];
    end
end

end