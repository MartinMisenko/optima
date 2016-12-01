% This MATLAB function represents non-convex optimization problem solver which
% implements Barrier and Interior Point methods with internal GDM
% algorithm. The function expects the optimization problem in general form:
% 
%       min  f(x)
%       s.t. h(x) <= 0
%  
% where f(x) and also h(x) could be convex or non-convex symbolic functions. 
% If we consider more constraints, the h(x) will be matrix (column vector)
% of symbolic expresions which represent left side of constraints.
% 
% Input arguments:
%   
%    fx    - (cost) function with symbolic variables
%    hx    - matrix (column vector) of symbolic expresions which represent 
%            left side of constraints
%    x0    - Iniatial point (in this phase FEASIBLE!!!)
%    kmax  - maximal No. of iteration for internal GDM
%    t     - coeficient for logarithmic approximation
%    tbt   - step size for GDM
%    alpha - GDM backtracking coeficient
%    beta  - GDM backtracking coeficient
%    eps   - precision
%    mi    - reducing factor of t
%    vp    - vector of used symbolic variables
%
%  © Simon Koniar -> sùaûnosti adresujte na: simon.koniar@gmail.com :)
function xx = IP_GDM_convex_feasible_start(fx, hx, kmax, t, tbt, alpha, beta,eps,mi,vp) 

%---------------Feasible starting point-------------------------------------
nx = size(vp,2);                        % number of variables
nc = size(hx,1);                        % number of constaints
con = matlabFunction(hx);               % convert constraint symbolic expression to an anonymous matlab function

s = sdpvar(1,1);                        % new variable
xs = sdpvar(nx,1);                      % initialization of yalmip variable
xs = num2cell(xs);                      % we need to separate previous yalmip variable into independent variables 

% min s
%  s.t. f_i(x)<=s                       % f_i(x) is i constraint

info = optimize(([con(xs{:})<=s]),s);

Feasible_start = [];                    % feasible starting point
for k = 1:nx
    Feasible_start = [Feasible_start;value(xs{k})]
end
%---------------Interior point - GDM---------------------------------------
xx = Feasible_start;
phi = 0;
for k = 1:nc
    phi = phi + log(hx(k));
end

while t > eps
    obj = [fx - t*phi];
    grad = [gradient(obj)]; 
    for k = 1:kmax 
        su = grad;
        for i = 1:nx
            su = subs(su,vp(i),xx(i)); % !!! tu je problem
        end
        dx = -real(double(su)); 
        if norm(dx)<eps, break, end 
        su1 = obj;
        su2 = obj;
        for i = 1:nx
            su1 = subs(su1,vp(i),xx(i)+tbt*dx(i));
            su2 = subs(su2,vp(i),xx(i));
        end
        a = su1;
        b = su2;
        while real(double(a))>real(double(b))+alpha*tbt*-dx'*dx 
            su1 = obj;
            su2 = obj;
            for i = 1:nx
                su1 = subs(su1,vp(i),xx(i)+tbt*dx(i));
                su2 = subs(su2,vp(i),xx(i));
            end
            a = su1;
            b = su2;
            tbt = beta*tbt;
        end 
        xx = xx + tbt*dx; 
    end
    xx =  double(xx)% po zmazani bodkociarky zobrazi vsetky priebezne vysledky
    t = t*mi;
end
end



