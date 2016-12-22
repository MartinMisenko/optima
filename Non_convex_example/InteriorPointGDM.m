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
%  © Simon Koniar -> sažnosti adresujte na: simon.koniar@gmail.com :)

function [xx,xy] = InteriorPointGDM(fx, hx, x0, kmax, t, tbt, alpha, beta,eps,mi,vp) 
xy = [];
nc = size(hx,1);
nx = size(x0,1);

% sem pride FAZA 0

xx = x0; % sem pride riesenie FAZY 0
phi = 0;
for k = 1:nc
    phi = phi + log(hx(k));
end

while t > eps
    obj = fx - t*phi;
    grad = gradient(obj); 
    for k = 1:kmax 
        dx = -real(double(subs(grad,vp',xx))); 
        if norm(dx)<eps, break, end 
        a = subs(obj,vp',xx+tbt*dx);
        b = subs(obj,vp',xx);
        while real(double(a))>real(double(b))+alpha*tbt*-dx'*dx 
            a = subs(obj,vp',xx+tbt*dx);
            b = subs(obj,vp',xx);
            tbt = beta*tbt;
        end 
        xx = xx + tbt*dx; 
    end
    xx =  double(xx)% po zmazani bodkociarky zobrazi vsetky priebezne vysledky
    xy = [xy,xx];
    t = t*mi;
end
 
end



