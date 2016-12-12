%% test

clear
close all

%% YALMIP problem

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);

obj = 100*(x2-x1^2)^2 + (1-x1)^2;
constr = [(x1-1/3)^2 + (x2-1/3)^2 - (1/3)^2 <= 0;... 
    0 <= x1 <= 0.5; 0.2 <= x2 <= 0.8];

sol = optimize(constr,obj);
value(x1),value(x2)

%% YALMIP phase 1

% sdpvar s1
% sdpvar s2
% sdpvar s3
% sdpvar s4
% sdpvar s5
s = sdpvar(5,1)

% obj = s1 + s2 + s3 + s4 + s5;
obj = sum(s);
% constr = [(x1-1/3)^2 + (x2-1/3)^2 - (1/3)^2 <= s1;...
%     -x1 <= s2;...
%     x1-0.5 <= s3;...
%     -x2+0.2<= s4;...
%     x2-0.8 <= s5;...
%     s1>= 0; s2>=0; s3>= 0; s4>=0; s5>=0];
constr = [(x1-1/3)^2 + (x2-1/3)^2 - (1/3)^2 <= s(1);...
    -x1 <= s(2);...
    x1-0.5 <= s(3);...
    -x2+0.2<= s(4);...
    x2-0.8 <= s(5);...
    s>=0];

sol = optimize(constr,obj);
% value(s1),value(s2),value(s3),value(s4),value(s5)
value(s)
value(x1),value(x2)


%% GDM problem

syms x11 x22 real

fx = 100*(x22-x11^2)^2 + (1-x11)^2;
hx = [(x11-1/3)^2 + (x22-1/3)^2 - (1/3)^2;...
    x11-0.5;...
    -x11;...
    x22-0.8;...
    -x22+0.2];
vp = [x11,x22];


%% PHASE 1

x0 = feasibleX(hx,vp)

%% GDM

t = 0.5;
kmax = 30;
tbt = 0.1;
alpha = 0.2;
beta = 0.5;
eps = 1e-6;
mi = 0.5;

xx = InteriorPointGDM(fx, hx, x0, kmax, t, tbt, alpha, beta,eps,mi,vp) 


