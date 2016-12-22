function x = feasibleX(hx,vp)

nx = length(vp); % number of variables
nc = size(hx,1); % number of constaints
con = matlabFunction(hx); % convert constraint symbolic expression to an anonymous matlab function

s = sdpvar(nc,1);                      
xs = sdpvar(nx,1);
xs = num2cell(xs); 

info = optimize(([con(xs{:})<=s; s>= 0]),sum(s));

x = [];                    % feasible starting point
for k = 1:nx
    x = [x;value(xs{k})];
end

end